# Ellenberg R and pH model
library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(patchwork)
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('stan_utility.R', local=util)
c_dark_trans <- "#80808080"
c_yellow_trans <- "#FFFF0080"

mang_cols <- unname(palette.colors()[c(2,6)])


# taking PH_long and X_Ell from 01b script
ph_ell_long <- full_join(X_Ell, PH_long) %>%
  filter(!Year %in% c(1990,2016)) %>%
  mutate(SQNUM = sapply(strsplit(REP_ID, "[A-Z]"), "[",1),
         YR = as.factor(Year)) %>%
  mutate(YRnum = as.integer(YR)) 
summary(ph_ell_long)

janitor::get_dupes(ph_ell_long, Year, REP_ID) # no duplicates

# check to see if there are patterns of missingness indicating issues with data
# combination
mice::md.pattern(select(ph_ell_long, Year, REP_ID, SM_R_UW, SM_R_W,
                        WH_R_UW, WH_R_W, pH, pH_CaCl2))

filter(ph_ell_long, is.na(pH)) %>%
  select(Year, REP_ID, SM_R_UW, SM_R_W,
         WH_R_UW, WH_R_W) %>%
  .$Year %>% table()

filter(ph_ell_long, Year == 1998) %>%
  select(Year, REP_ID, SM_R_UW, SM_R_W,
         WH_R_UW, WH_R_W, pH) %>%
  summary()
# < half X plots in 1998 have matching soil pH data

# Linear model ####
mod_data <- ph_ell_long %>% ungroup() %>%
  mutate(Ell_R = as.numeric(scale(WH_R_W)),
         Soil_pH = as.numeric(scale(pH))) %>%
  select(Soil_pH, Ell_R, Year, YR, YRnum, SQNUM, REP_ID) %>%
  filter(!(is.na(Soil_pH) | is.na(Ell_R)))
get_prior(Ell_R ~ Soil_pH*YR + (1|SQNUM/REP_ID), data = mod_data,
          family = gaussian())

mod_pr <- c(prior(normal(0,1), class = "b"),
            prior(student_t(3, 0, 1), class = "Intercept"),
            prior(student_t(3, 0, 1), class = "sd"),
            prior(student_t(3, 0, 1), class = "sigma"))

mod_ln <- brm(Ell_R ~ Soil_pH*YR + (1|SQNUM/REP_ID), 
              data = mod_data,
              prior = mod_pr,
              cores = 4)
summary(mod_ln)
plot(mod_ln)
plot(conditional_effects(mod_ln), points = TRUE)

# autoregressive model
get_prior(Ell_R ~ Soil_pH + (1|SQNUM) +
            ar(time = YRnum, gr = REP_ID), data = mod_data)

mod_pr <- c(prior(normal(0.3,0.5), class = "b"),
            prior(normal(0.2,0.1), class = "ar"),
            prior(student_t(3, 0, 1), class = "Intercept"),
            prior(student_t(3, 0, 1), class = "sd"),
            prior(student_t(3, 0, 1), class = "sigma"))

mod_ln <- brm(Ell_R ~ Soil_pH + (1|SQNUM) +
                ar(time = YRnum, gr = REP_ID), 
              data = mod_data, iter = 4000,
              cores = 4, prior = mod_pr,
              file = "Outputs/Models/EllR_pH_modln_020920")
saveRDS(mod_ln, "Outputs/Models/EllR_pH_modln_020920.rds")
summary(mod_ln)
plot(mod_ln)
plot(conditional_effects(mod_ln), points = TRUE)


# Non-linear model linking pH and Ellenberg R ####

# Assymetrical sigmoidal curve, with impact of pH varying by year
# c1 is the max y (minus c5)
# c2 is the rate of increase for the curve
# c3 is the x at curve midpoint
# c4 is the amount of asymmetry
# c5 is the min y
# c6 is the random effect(square shifts whole curve up/down)
mod_data$YR <- relevel(mod_data$YR, "2007")
get_prior(bf(Ell_R ~ c1/((1 + exp(-c2*(Soil_pH - c3)))^c4) + c5, 
             c5 ~ (1|SQNUM), c1 ~ 1, c2 + c3 + c4 ~ YR, nl= TRUE),
          data = mod_data)

pr <- c(prior(normal(2,0.2), nlpar = "c1"),
        prior(normal(1.5,0.2), nlpar = "c2"),
        prior(normal(0,0.1), nlpar = "c2", coef = "YR1978"),
        prior(normal(0,0.1), nlpar = "c2", coef = "YR1998"),
        prior(normal(0,0.1), nlpar = "c2", coef = "YR2019"),
        prior(normal(-0.5,1), nlpar = "c3", coef = "Intercept"),
        prior(normal(0,0.1), nlpar = "c3", coef = "YR1978"),
        prior(normal(0,0.1), nlpar = "c3", coef = "YR1998"),
        prior(normal(0,0.1), nlpar = "c3", coef = "YR2019"),
        prior(gamma(2,2), nlpar = "c4", lb = 0),
        prior(normal(0,0.1), nlpar = "c4", coef = "YR1978"),
        prior(normal(0,0.1), nlpar = "c4", coef = "YR1998"),
        prior(normal(0,0.1), nlpar = "c4", coef = "YR2019"),
        prior(normal(-1.5,0.2), nlpar = "c5", coef = "Intercept"),
        prior(student_t(3,0,.1), nlpar = "c5", group = "SQNUM", class = "sd"))

# prior simulation
mod_pr_only <- brm(bf(Ell_R ~ c1/((1 + exp(-c2*(Soil_pH - c3)))^c4) + c5, 
                      c5 ~ (1|SQNUM), c1 ~ 1, c2 + c3 + c4 ~ YR, nl= TRUE),
                   data = mod_data, prior = pr, cores = 4, chains = 4,
                   sample_prior = "only")
summary(mod_pr_only)
plot(mod_pr_only)
plot(conditional_effects(mod_pr_only))

pp_check(mod_pr_only)

pred_modpr <- predict(mod_pr_only)

plot(mod_data$Ell_R, pred_modpr[,1])


# with data
mod_nl <- brm(bf(Ell_R ~ c1/((1 + exp(-c2*(Soil_pH - c3)))^c4) + c5, 
                 c5 ~ (1|SQNUM), c1 ~ 1, c2 + c3 + c4 ~ YR, nl= TRUE),
              data = mod_data, prior = pr, cores = 4, chains = 4)
summary(mod_nl)
plot(mod_nl)
plot(conditional_effects(mod_pr_only))

pp_check(mod_pr_only)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# DIFFERENCE MODELS ####
# Ellenberg and pH difference models

# check which deposition statistic better predicts pH and soil N
phdep_data <- PH %>%
  rename_with(~gsub("diff","PH_diff",.x)) %>%
  mutate(PHC_diff0719 = PHC_2019 - PHC_2007) %>%
  select(REP_ID, contains("diff7898"), contains("diff9807"),
         contains("diff0719")) %>%
  pivot_longer(contains("diff"), names_to = c("pH","Time_period"),
               names_sep= "_diff") %>%
  na.omit()  %>%
  pivot_wider(id_cols = c(REP_ID, Time_period), names_from = pH) %>%
  mutate(Year = recode(Time_period,
                       "7898" = 1998,
                       "9807" = 2007,
                       "0719" = 2019)) %>%
  full_join(CN) %>%
  left_join(select(CS_plot_atdep, -Habitat)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Sdep = as.numeric(scale(Sdep)),
         Ndep = as.numeric(scale(Ndep)))

phdep_8yr <- filter(phdep_data,
                    Duration == "8yr")
mod_pr <- c(prior(normal(0,1), class = "b"),
            prior(normal(0,0.5), class = "Intercept"),
            prior(normal(0,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))
mod_ph_sdep8yr <- brm(PH ~ Sdep + (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID),
                      family = "student", prior = mod_pr, data = phdep_8yr,
                      cores = 4, iter = 4000)
summary(mod_ph_sdep8yr)

phdep_20yr <- filter(phdep_data,
                     Duration == "20yr")
mod_ph_sdep20yr <- brm(PH ~ Sdep + (1|SQUARE) +
                         ar(time = YRnm, gr = REP_ID),
                       family = "student", prior = mod_pr, data = phdep_20yr,
                       cores = 4, iter = 4000)
summary(mod_ph_sdep20yr)
plot(conditional_effects(mod_ph_sdep20yr), points = TRUE)
plot(conditional_effects(mod_ph_sdep8yr), points = TRUE)

sdep_eff <- 
  rbind(cbind(as.data.frame(fixef(mod_ph_sdep8yr)),
              Duration = c("8yr","8yr")),
        cbind(as.data.frame(fixef(mod_ph_sdep20yr)),
              Duration = c("20yr","20yr")))[c(2,4),]
ggplot(sdep_eff, aes(x = Duration)) +
  geom_pointrange(aes(y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  scale_y_continuous(limits = c(-0.5,0)) +
  theme(panel.grid.major = element_line(colour = "gray90")) +
  labs(x = "S deposition difference duration",
       y = "Estimated impact of S on pH")
ggsave("S dep duration and impact on pH.png",
       path= "Outputs/Models/", width = 12, height = 12,
       units = "cm")

# n deposition effect on NC
mod_pr <- c(prior(normal(0,1), class = "b"),
            prior(normal(0,0.5), class = "Intercept"),
            prior(normal(0,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))
mod_nc_ndep8yr <- brm(NC_RATIO ~ Ndep + (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID),
                      prior = mod_pr, data = phdep_8yr,
                      cores = 4, iter = 4000)
summary(mod_nc_ndep8yr)

mod_nc_ndep20yr <- update(mod_nc_ndep8yr, newdata = phdep_20yr,
                          cores = 4, iter = 4000)
summary(mod_nc_ndep20yr)
plot(conditional_effects(mod_nc_ndep20yr), points = TRUE)
plot(conditional_effects(mod_nc_ndep8yr), points = TRUE)

ndep_eff <- 
  rbind(cbind(as.data.frame(fixef(mod_nc_ndep8yr)),
              Duration = c("8yr","8yr")),
        cbind(as.data.frame(fixef(mod_nc_ndep20yr)),
              Duration = c("20yr","20yr")))[c(2,4),]
ggplot(ndep_eff, aes(x = Duration)) +
  geom_pointrange(aes(y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  scale_y_continuous(limits = c(0,0.016)) +
  theme(panel.grid.major = element_line(colour = "gray90")) +
  labs(x = "N deposition difference duration",
       y = "Estimated impact of N on N:C")
ggsave("N dep duration and impact on NC.png",
       path= "Outputs/Models/", width = 12, height = 12,
       units = "cm")

# N dep effect on pH
mod_ph_ndep8yr <- brm(PH ~ Ndep + (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID),
                      family = "student", prior = mod_pr, data = phdep_8yr,
                      cores = 4, iter = 4000)
summary(mod_ph_ndep8yr)

mod_ph_ndep20yr <- brm(PH ~ Ndep + (1|SQUARE) +
                         ar(time = YRnm, gr = REP_ID),
                       family = "student", prior = mod_pr, data = phdep_20yr,
                       cores = 4, iter = 4000)
summary(mod_ph_ndep20yr)

ndep_pheff <- 
  rbind(cbind(as.data.frame(fixef(mod_ph_ndep8yr)),
              Duration = c("8yr","8yr")),
        cbind(as.data.frame(fixef(mod_ph_ndep20yr)),
              Duration = c("20yr","20yr")))[c(2,4),]
ggplot(ndep_pheff, aes(x = Duration)) +
  geom_pointrange(aes(y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  scale_y_continuous(limits = c(0,0.4)) +
  theme(panel.grid.major = element_line(colour = "gray90")) +
  labs(x = "N deposition difference duration",
       y = "Estimated impact of N on pH")
ggsave("N dep duration and impact on pH.png",
       path= "Outputs/Models/", width = 12, height = 12,
       units = "cm")


# data manipulation check if just limiting to the 78-98, 98-07 and 07-19 periods
# is okay in terms of data coverage
Ell_R <- X_Ell %>%
  select(Year, REP_ID, contains("_R_")) %>%
  pivot_longer(contains("_R_"), names_to = "Ellenberg") %>%
  pivot_wider(id_cols = REP_ID, names_from = c(Ellenberg,Year)) %>%
  mutate(WH_R_W_diff7898 = WH_R_W_1998 - WH_R_W_1978,
         WH_R_W_diff9807 = WH_R_W_2007 - WH_R_W_1998,
         WH_R_W_diff0719 = WH_R_W_2019 - WH_R_W_2007,
         WH_R_W_diff7807 = WH_R_W_2007 - WH_R_W_1978,
         WH_R_W_diff7819 = WH_R_W_2019 - WH_R_W_1978,
         WH_R_W_diff9819 = WH_R_W_2019 - WH_R_W_1998,
  )
mice::md.pattern(select(Ell_R, contains("diff")))
# 33 missing 1998 but having 2007/1978, 4 missing 2007 but 
# having 2019 and 1998 - this seems fine to me to just use the 78/98, 
# 98/07 and 07/19 changes

# get Ellenberg R difference data 
Ell_R <- X_Ell %>%
  select(Year, REP_ID, contains("_R_")) %>%
  pivot_longer(contains("_R_"), names_to = "Ellenberg") %>%
  pivot_wider(id_cols = REP_ID, names_from = c(Ellenberg,Year)) %>%
  mutate(WH_R_W_diff7898 = WH_R_W_1998 - WH_R_W_1978,
         WH_R_W_diff9807 = WH_R_W_2007 - WH_R_W_1998,
         WH_R_W_diff0719 = WH_R_W_2019 - WH_R_W_2007,
         WH_R_UW_diff7898 = WH_R_UW_1998 - WH_R_UW_1978,
         WH_R_UW_diff9807 = WH_R_UW_2007 - WH_R_UW_1998,
         WH_R_UW_diff0719 = WH_R_UW_2019 - WH_R_UW_2007,
         SM_R_W_diff7898 = SM_R_W_1998 - SM_R_W_1978,
         SM_R_W_diff9807 = SM_R_W_2007 - SM_R_W_1998,
         SM_R_W_diff0719 = SM_R_W_2019 - SM_R_W_2007,
         SM_R_UW_diff7898 = SM_R_UW_1998 - SM_R_UW_1978,
         SM_R_UW_diff9807 = SM_R_UW_2007 - SM_R_UW_1998,
         SM_R_UW_diff0719 = SM_R_UW_2019 - SM_R_UW_2007
  ) %>%
  select(REP_ID, contains("diff")) %>%
  pivot_longer(contains("diff"), names_sep = "_diff",
               names_to = c("Ell","Time_period")) %>%
  na.omit() %>%
  pivot_wider(id_cols = c(REP_ID, Time_period), names_from = Ell) 

# combine Ellenberg R difference data with pH and climatic data
ELL_pH <- PH %>%
  rename_with(~gsub("diff","PH_diff",.x)) %>%
  mutate(PHC_diff0719 = PHC_2019 - PHC_2007) %>%
  select(REP_ID, contains("diff7898"), contains("diff9807"),
         contains("diff0719")) %>%
  pivot_longer(contains("diff"), names_to = c("pH","Time_period"),
               names_sep= "_diff") %>%
  na.omit()  %>%
  pivot_wider(id_cols = c(REP_ID, Time_period), names_from = pH) %>%
  full_join(Ell_R) %>%
  mutate(Year = recode(Time_period,
                       "7898" = 1998,
                       "9807" = 2007,
                       "0719" = 2019),
         Year1 = recode(Time_period,
                        "7898" = 1978,
                        "9807" = 1998,
                        "0719" = 2007)) %>%
  left_join(BH_IMP) %>%
  left_join(filter(CS_plot_atdep, Duration == "8yr") %>%
              select(REP_ID, Year, Ndep, Sdep) %>%
              mutate(Year = ifelse(Year == 2018, 2019, Year))) %>%
  left_join(cs_rainfall_diff) %>%
  left_join(cs_rainfall_averages) %>%
  left_join(ELL_QA_diff) %>%
  left_join(PH_QA_diff) %>%
  left_join(rename(PH_long, Year1 = Year, Year1_pH = pH, 
                   Year1_pHCaCl2 = pH_CaCl2)) %>%
  left_join(rename(PH_long, Year2_pH = pH, 
                   Year2_pHCaCl2 = pH_CaCl2)) %>%
  left_join(CN)

janitor::get_dupes(ELL_pH, Time_period, REP_ID) # no duplicates
mice::md.pattern(ELL_pH)
View(filter(ELL_pH, is.na(WH_R_UW) & !is.na(PH)))



# simple difference model
# model each Ellenberg R change as a function of pH change 
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W,
         Management = ifelse(Management == "High",1,0)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, Management) %>%
  na.omit()

# normal model response - to show not great
get_prior(Ell ~ PH*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
          data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(normal(0.4,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))
normal_mod <- brm(Ell ~ PH*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                  data = mod_data, prior = mod_pr, sample_prior = "only",
                  cores = 4, iter = 5000)
pp_check(normal_mod, nsamples = 50) + scale_x_continuous(limits =c(-5,5))

# prior checks
get_prior(Ell ~ PH*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
          family = "student", data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(gamma(2,0.1), class = "nu"),
            prior(normal(0.4,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))

# prior predictive check
prior_mod <- brm(Ell ~ PH*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                 family = "student", sample_prior = "only",
                 data = mod_data, prior = mod_pr, cores = 4, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50)
# underestimates peak at 0 on average but gets close

# ellenberg R model - weighted Ellenberg R for whole plot
ell_ph_diff <- brm(Ell ~ PH*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                   data = mod_data, prior = mod_pr, family = "student",
                   cores = 4, iter = 5000, control = list(adapt_delta = 0.99),
                   file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_WHW_PH_HAB")
summary(ell_ph_diff)
plot(ell_ph_diff, ask = FALSE)
pp_check(ell_ph_diff)

# unweighted ellenberg r whole plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW,
         Management = ifelse(Management == "High",1,0)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, Management) %>%
  na.omit() 
ell_ph_diffuw <- update(ell_ph_diff, newdata = mod_data, cores=4, iter = 5000,
                        file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_WHUW_PH_HAB",
                        control = list(adapt_delta = 0.95))
summary(ell_ph_diffuw)
plot(ell_ph_diffuw, ask = FALSE)
pp_check(ell_ph_diffuw)


# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W,
         Management = ifelse(Management == "High",1,0)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, Management) %>%
  na.omit() 
ell_ph_diffsmw <- update(ell_ph_diff, newdata = mod_data, cores = 4, iter = 5000,
                         file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_SMW_PH_HAB",
                         control = list(adapt_delta = 0.95))
summary(ell_ph_diffsmw)
plot(ell_ph_diffsmw, ask = FALSE)
pp_check(ell_ph_diffsmw)

# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW,
         Management = ifelse(Management == "High",1,0)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, Management) %>%
  na.omit() 
ell_ph_diffsmuw <- update(ell_ph_diff, newdata = mod_data, cores = 4, iter = 5000,
                          control = list(adapt_delta = 0.95),
                          file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_SMUW_PH_HAB")
summary(ell_ph_diffsmuw)
plot(ell_ph_diffsmuw, ask = FALSE)
pp_check(ell_ph_diffsmuw)


# Plot for comparing pH effects on different Ellenberg R scores
nd <- 
  tibble(PH = seq(from = -3.25, to = 3.5, length.out = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 30),
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(ell_ph_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(ell_ph_diffuw, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(ell_ph_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(ell_ph_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "pH", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("pH effect on Ellenberg R score versions.png",
       path = "Outputs/Models/Difference/Univariate_nomeaserror", 
       width = 16, height = 12, units = "cm")

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = PH, y = value, colour = Management),
             alpha = 0.2) +
  geom_smooth(data = f,
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "pH", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("pH effect on Ellenberg R score versions with data.png",
       path = "Outputs/Models/Difference/Univariate_nomeaserror", 
       width = 16, height = 12, units = "cm")


# ** pH diw vs cacl2 comparison ####
# only one time transition so no temporal autocorrelation included
mod_data <- ELL_pH %>%
  mutate(SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W,
         Management = ifelse(Management == "High",1,0)) %>%
  select(REP_ID, SQUARE, Ell, PH, PHC, Management) %>%
  na.omit()

get_prior(Ell ~ PH*Management + (1|SQUARE),
          data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))

ell_ph_diff <- brm(Ell ~ PH*Management + (1|SQUARE), family = "student",
                   data = mod_data, prior = mod_pr, cores = 4, iter = 5000, 
                   file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_WHW_PH_HAB_0719only",
                   save_pars = save_pars(all = TRUE), control = list(adapt_delta = 0.95))
summary(ell_ph_diff)
plot(ell_ph_diff, ask = FALSE)
pp_check(ell_ph_diff)
ell_ph_diff <- add_criterion(ell_ph_diff, "loo", moment_match = TRUE, reloo = TRUE,
                             file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_WHW_PH_HAB_0719only")

ell_phc_diff <- brm(Ell ~ PHC*Management + (1|SQUARE), family = "student",
                    data = mod_data, prior = mod_pr, cores = 4, iter = 5000, 
                    file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_WHW_PHC_HAB_0719only",
                    save_pars = save_pars(all = TRUE),  control = list(adapt_delta = 0.95))
summary(ell_phc_diff)
plot(ell_phc_diff, ask = FALSE)
pp_check(ell_phc_diff)
ell_phc_diff <- add_criterion(ell_phc_diff, "loo", moment_match = TRUE, reloo = TRUE,
                              file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_WHW_PHC_HAB_0719only")

loo_compare(ell_ph_diff, ell_phc_diff)
#               elpd_diff se_diff
# ell_phc_diff  0.0       0.0   
# ell_ph_diff  -0.5       1.5

# unweighted ellenberg r whole plot
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW,
         Management = ifelse(Management == "High",1,0)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, PHC, Management) %>%
  na.omit()
ell_ph_diffuw <- update(ell_ph_diff, newdata = mod_data, cores = 4, iter = 5000,
                        file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_WHUW_PH_HAB_0719only")
summary(ell_ph_diffuw)
plot(ell_ph_diffuw, ask = FALSE)
pp_check(ell_ph_diffuw)
ell_ph_diffuw <- add_criterion(ell_ph_diffuw, "loo", moment_match = TRUE, reloo = TRUE,
                               file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_WHUW_PH_HAB_0719only")

ell_phc_diffuw <- update(ell_phc_diff, newdata = mod_data, cores = 4, iter = 5000,
                         file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_WHUW_PHC_HAB_0719only")
summary(ell_phc_diffuw)
plot(ell_phc_diffuw, ask = FALSE)
pp_check(ell_phc_diffuw)
ell_phc_diffuw <- add_criterion(ell_phc_diffuw, "loo", moment_match = TRUE, reloo = TRUE,
                                file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_WHUW_PHC_HAB_0719only")

loo_compare(ell_ph_diffuw, ell_phc_diffuw)
# elpd_diff se_diff
# ell_phc_diffuw  0.0       0.0   
# ell_ph_diffuw  -0.2       1.1  

# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W,
         Management = ifelse(Management == "High",1,0)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, PHC, Management) %>%
  na.omit()
ell_ph_diffsmw <- update(ell_ph_diff, newdata = mod_data, cores = 4, iter = 5000,
                         file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_SMW_PH_HAB_0719only")
summary(ell_ph_diffsmw)
plot(ell_ph_diffsmw, ask = FALSE)
pp_check(ell_ph_diffsmw)
ell_ph_diffsmw <- add_criterion(ell_ph_diffsmw, "loo", moment_match = TRUE, reloo = TRUE,
                                file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_SMW_PH_HAB_0719only")

ell_phc_diffsmw <- update(ell_phc_diff, newdata = mod_data, cores = 4, iter = 5000,
                          file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_SMW_PHC_HAB_0719only")
summary(ell_phc_diffsmw)
plot(ell_phc_diffsmw, ask = FALSE)
pp_check(ell_phc_diffsmw)
ell_phc_diffsmw <- add_criterion(ell_phc_diffsmw, "loo", moment_match = TRUE, reloo = TRUE,
                                 file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_SMW_PHC_HAB_0719only")

loo_compare(ell_ph_diffsmw, ell_phc_diffsmw)
# elpd_diff se_diff
# ell_ph_diffsmw   0.0       0.0   
# ell_phc_diffsmw -0.7       1.9   

# unweighted ellenberg r small plot
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW,
         Management = ifelse(Management == "High",1,0)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, PHC, Management) %>%
  na.omit()
ell_ph_diffsmuw <- update(ell_ph_diff, newdata = mod_data, cores = 4, iter = 5000,
                          file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_SMUW_PH_HAB_0719only")
summary(ell_ph_diffsmuw)
plot(ell_ph_diffsmuw, ask = FALSE)
pp_check(ell_ph_diffsmuw)
ell_ph_diffsmuw <- add_criterion(ell_ph_diffsmuw, "loo", moment_match = TRUE, reloo = TRUE,
                                 file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_SMUW_PH_HAB_0719only")

ell_phc_diffsmuw <- update(ell_phc_diff, newdata = mod_data, cores = 4, iter = 5000,
                           file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_SMUW_PHC_HAB_0719only")
summary(ell_phc_diffsmuw)
plot(ell_phc_diffsmuw, ask = FALSE)
pp_check(ell_phc_diffsmuw)
ell_phc_diffsmuw <- add_criterion(ell_phc_diffsmuw, "loo", moment_match = TRUE, reloo = TRUE,
                                  file = "Outputs/Models/Difference/pH DIW vs pH CaCl2/EllR_SMUW_PHC_HAB_0719only")

loo_compare(ell_ph_diffsmuw, ell_phc_diffsmuw)
# elpd_diff se_diff
# ell_ph_diffsmuw  0.0       0.0    
# ell_phc_diffsmuw 0.0       1.1

# Change in pH (DIW) and change in pH (CaCl2) equivalent predictors of Ellenberg
# R change. 

bayes_R2(ell_ph_diff) #0.087
bayes_R2(ell_ph_diffuw) #0.174
bayes_R2(ell_ph_diffsmw) #0.098
bayes_R2(ell_ph_diffsmuw) #0.047

bayes_R2(ell_phc_diff) #0.088
bayes_R2(ell_phc_diffuw) #0.181
bayes_R2(ell_phc_diffsmw) #0.099
bayes_R2(ell_phc_diffsmuw) #0.045

# Strength of relationship for both pH DIW and CaCl2:
# large unweighted > small weighted > large weighted > small unweighted

# Atmospheric Deposition ####
## ** Sdep ####
# model each Ellenberg R change as a function of Sdep
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W,
         Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Sdep, Management) %>%
  na.omit()

# prior checks
get_prior(Ell ~ Sdep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
          family = "student", data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(gamma(2,0.1), class = "nu"),
            prior(normal(0.4,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))

# prior predictive check model
prior_mod <- brm(Ell ~ Sdep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                 family = "student", sample_prior = "only",
                 data = mod_data, prior = mod_pr, cores = 4, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50)
# underestimates peak at 0 on average but gets close

# ellenberg R model - weighted Ellenberg R for whole plot
ell_Sdep_diff <- brm(Ell ~ Sdep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                     data = mod_data, prior = mod_pr, family = "student",
                     cores = 4, iter = 5000, 
                     file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_WHW_Sdep_HAB")
summary(ell_Sdep_diff)
plot(ell_Sdep_diff, ask = FALSE)
pp_check(ell_Sdep_diff)

# unweighted ellenberg r whole plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW,
         Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Sdep, Management) %>%
  na.omit() 
ell_Sdep_diffuw <- update(ell_Sdep_diff, newdata = mod_data, cores=4, iter = 5000,
                          file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_WHUW_Sdep_HAB")
summary(ell_Sdep_diffuw)
plot(ell_Sdep_diffuw, ask = FALSE)
pp_check(ell_Sdep_diffuw)


# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W,
         Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Sdep, Management) %>%
  na.omit() 
ell_Sdep_diffsmw <- update(ell_Sdep_diff, newdata = mod_data, cores=4, iter = 5000,
                           control = list(adapt_delta = 0.99),
                           file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_SMW_Sdep_HAB")
summary(ell_Sdep_diffsmw)
plot(ell_Sdep_diffsmw, ask = FALSE)
pp_check(ell_Sdep_diffsmw)

# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW,
         Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Sdep, Management) %>%
  na.omit() 
ell_Sdep_diffsmuw <- update(ell_Sdep_diff, newdata = mod_data, cores=4, iter = 5000,
                            control = list(adapt_delta = 0.99),
                            file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_SMUW_Sdep_HAB")
summary(ell_Sdep_diffsmuw)
plot(ell_Sdep_diffsmuw, ask = FALSE)
pp_check(ell_Sdep_diffsmuw)


# Plot for comparing Sdep effects on different Ellenberg R scores
nd <- 
  tibble(Sdep = seq(from = -5.5, to = 2, length.out = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 30),
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(ell_Sdep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(ell_Sdep_diffuw, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(ell_Sdep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(ell_Sdep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W, WH_R_UW, SM_R_W, SM_R_UW,
         Sdep, Management, REP_ID, Time_period) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Sdep = as.numeric(scale(Sdep)),
         Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "S deposition", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("Sdep effect on Ellenberg R score versions.png",
       path = "Outputs/Models/Difference/Univariate_nomeaserror", 
       width = 16, height = 12, units = "cm")

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = Sdep, y = value, colour = Management),
             alpha = 0.2) +
  geom_smooth(data = f,
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Sdep", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("Sdep effect on Ellenberg R score versions with data.png",
       path = "Outputs/Models/Difference/Univariate_nomeaserror", 
       width = 16, height = 12, units = "cm")

# ** Ndep ####
# model each Ellenberg R change as a function of Ndep
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W,
         Management = ifelse(Management == "High",1,0),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ndep, Management) %>%
  na.omit()

# prior checks
get_prior(Ell ~ Ndep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
          family = "student", data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(gamma(2,0.1), class = "nu"),
            prior(normal(0.4,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))

# prior predictive check model
prior_mod <- brm(Ell ~ Ndep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                 family = "student", sample_prior = "only",
                 data = mod_data, prior = mod_pr, cores = 4, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50)
# underestimates peak at 0 on average but gets close

# ellenberg R model - weighted Ellenberg R for whole plot
ell_Ndep_diff <- brm(Ell ~ Ndep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                     data = mod_data, prior = mod_pr, family = "student",
                     cores = 4, iter = 5000, 
                     file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_WHW_Ndep_HAB")
summary(ell_Ndep_diff)
plot(ell_Ndep_diff, ask = FALSE)
pp_check(ell_Ndep_diff)

# unweighted ellenberg r whole plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW,
         Management = ifelse(Management == "High",1,0),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ndep, Management) %>%
  na.omit() 
ell_Ndep_diffuw <- update(ell_Ndep_diff, newdata = mod_data, cores=4, iter = 5000,
                          file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_WHUW_Ndep_HAB")
summary(ell_Ndep_diffuw)
plot(ell_Ndep_diffuw, ask = FALSE)
pp_check(ell_Ndep_diffuw)


# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W,
         Management = ifelse(Management == "High",1,0),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ndep, Management) %>%
  na.omit() 
ell_Ndep_diffsmw <- update(ell_Ndep_diff, newdata = mod_data, cores = 4, iter = 5000,
                           control = list(adapt_delta = 0.99),
                           file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_SMW_Ndep_HAB")
summary(ell_Ndep_diffsmw)
plot(ell_Ndep_diffsmw, ask = FALSE)
pp_check(ell_Ndep_diffsmw)

# unweighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW,
         Management = ifelse(Management == "High",1,0),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ndep, Management) %>%
  na.omit() 
ell_Ndep_diffsmuw <- update(ell_Ndep_diff, newdata = mod_data, cores = 4, iter = 5000,
                            control = list(adapt_delta = 0.99),
                            file = "Outputs/Models/Difference/Univariate_nomeaserror/EllR_SMUW_Ndep_HAB")
summary(ell_Ndep_diffsmuw)
plot(ell_Ndep_diffsmuw, ask = FALSE)
pp_check(ell_Ndep_diffsmuw)


# Plot for comparing Ndep effects on different Ellenberg R scores
nd <- 
  tibble(Ndep = seq(from = -2, to = 5.5, length.out = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 30),
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(ell_Ndep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(ell_Ndep_diffuw, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(ell_Ndep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(ell_Ndep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W, WH_R_UW, SM_R_W, SM_R_UW,
         Ndep, Management, REP_ID, Time_period) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Ndep = as.numeric(scale(Ndep)),
         Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "N deposition", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("Ndep effect on Ellenberg R score versions.png",
       path = "Outputs/Models/Difference/Univariate_nomeaserror", 
       width = 16, height = 12, units = "cm")

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = Ndep, y = value, colour = Management),
             alpha = 0.2) +
  geom_smooth(data = f,
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Ndep", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("Ndep effect on Ellenberg R score versions with data.png",
       path = "Outputs/Models/Difference/Univariate_nomeaserror", 
       width = 16, height = 12, units = "cm")


# pH
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep)),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, PH, Sdep, Ndep, Management) %>%
  na.omit() 
ph_Sdep_diff <- brm(PH ~ Sdep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                    data = mod_data, prior = mod_pr, family = "student",
                    cores = 4, iter = 5000, 
                    file = "Outputs/Models/Difference/Univariate_nomeaserror/PH_Sdep_HAB")
summary(ph_Sdep_diff)
plot(ph_Sdep_diff, ask = FALSE)
pp_check(ph_Sdep_diff)

ph_Ndep_diff <- brm(PH ~ Ndep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                    data = mod_data, prior = mod_pr, family = "student",
                    cores = 4, iter = 5000, 
                    file = "Outputs/Models/Difference/Univariate_nomeaserror/PH_Ndep_HAB")
summary(ph_Ndep_diff)
plot(ph_Ndep_diff, ask = FALSE)
pp_check(ph_Ndep_diff)


# * Measurement error ####

# The following code (up until the graphs) was run in parallel on a modelling PC with
# many cores
library(brms)
library(dplyr)
options(future.globals.maxSize = 10e8)
library(future)
library(future.apply)
plan(multisession)

ELL_pH <- read.csv("Outputs/ELL_pH_181120.csv")

mod_pr <- mod_pr <- c(prior(normal(0,0.5), class = "b"),
                      prior(normal(0,0.25), class = "Intercept"),
                      prior(normal(2,0.5), class = "ar"),
                      prior(student_t(3,0,1), class = "sd"),
                      prior(student_t(3,0,1), class = "sigma"))

mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Response = PH) %>%
  select(REP_ID, SQUARE, YRnm, Management, Sdep, Ndep,
         Response,
         RESP_SE = PH_DIW_SE_NORM) %>%
  mutate(Predictor = as.numeric(scale(Sdep))) %>%
  na.omit()

set.seed(78342835)
init_mod <- brm(Response | mi(RESP_SE) ~ Predictor*Management +
                  (1|p|SQUARE) +
                  ar(time = YRnm, gr = REP_ID),
                data = mod_data, prior = mod_pr,
                save_pars = save_pars(all = TRUE, latent = TRUE), 
                file = "PH__Sdep_HAB",
                cores = 4, iter = 20000, thin = 4,
                seed = 78342835)

data_to_run <- data.frame(
  Response = rep(c("WH_R_W","WH_R_UW","SM_R_W","SM_R_UW","PH"),
                 each = 2),
  RESP_SE = rep(c("ELL_WH_W_SE_NORM","ELL_WH_UW_SE_NORM",
                  "ELL_SM_W_SE_NORM","ELL_SM_UW_SE_NORM",
                  "PH_DIW_SE_NORM"),
                each = 2),
  Predictor = rep(c("Sdep","Ndep"),5)
)

mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1))

# Function for updating the initial model with the new response/predictor combo
mod_update <- function(Response, RESP_SE, Predictor, mod_data,init_mod){
  mod_data$Response <- mod_data[,Response]
  mod_data$RESP_SE <- mod_data[,RESP_SE]
  mod_data$Predictor <- mod_data[,Predictor]
  
  mod_data <- select(mod_data, Response, Predictor, RESP_SE,
                     Management, YRnm, SQUARE, REP_ID)
  
  filename <- paste(Response, Predictor, sep = "__")
  
  fit <- update(init_mod, newdata = mod_data,
                cores = 4, iter = 20000, thin = 4,
                seed = 78342835, file = filename)
}

# futures code for running this on many cores
set.seed(78342835)
testrun <- future_mapply(mod_update, 
                         data_to_run$Response,
                         data_to_run$RESP_SE,
                         data_to_run$Predictor,
                         list(mod_data, mod_data,mod_data,mod_data,mod_data,
                              mod_data, mod_data,mod_data,mod_data,mod_data), 
                         list(init_mod,init_mod,init_mod,init_mod,init_mod,
                              init_mod,init_mod,init_mod,init_mod,init_mod),
                         future.seed = TRUE, SIMPLIFY = FALSE)

# All of the Ellenberg ~ Sdep models had divergent transitions
pars_plot <- c("b_Intercept","b_Predictor","b_ManagementLow",
               "b_Predictor:ManagementLow","ar[1]","sd_SQUARE__Intercept",
               "sigma","Intercept" )
pairs(testrun[[3]]$fit, pars = pars_plot, log = TRUE)
pairs(testrun[[5]]$fit, pars = pars_plot, log = TRUE)
pairs(testrun[[7]]$fit, pars = pars_plot, log = TRUE)

# Plot for comparing Sdep effects on different Ellenberg R scores
nd <- 
  tibble(Sdep = seq(from = -5.5, to = 2, length.out = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 30),
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(ell_Sdep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(ell_Sdep_diffuw, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(ell_Sdep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(ell_Sdep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W, WH_R_UW, SM_R_W, SM_R_UW,
         Sdep, Management, REP_ID, Time_period) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Sdep = as.numeric(scale(Sdep)),
         Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "S deposition", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("Sdep effect on Ellenberg R score versions.png",
       path = "Outputs/Models/Difference/Univariate_measerror", 
       width = 16, height = 12, units = "cm")

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = Sdep, y = value, colour = Management),
             alpha = 0.2) +
  geom_smooth(data = f,
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Sdep", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("Sdep effect on Ellenberg R score versions with data.png",
       path = "Outputs/Models/Difference/Univariate_measerror", 
       width = 16, height = 12, units = "cm")



# Plot for comparing Ndep effects on different Ellenberg R scores
nd <- 
  tibble(Ndep = seq(from = -2, to = 5.5, length.out = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 30),
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(ell_Ndep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(ell_Ndep_diffuw, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(ell_Ndep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(ell_Ndep_diff, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W, WH_R_UW, SM_R_W, SM_R_UW,
         Ndep, Management, REP_ID, Time_period) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Ndep = as.numeric(scale(Ndep)),
         Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "N deposition", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("Ndep effect on Ellenberg R score versions.png",
       path = "Outputs/Models/Difference/Univariate_measerror", 
       width = 16, height = 12, units = "cm")

ggplot(plot_dat) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = Ndep, y = value, colour = Management),
             alpha = 0.2) +
  geom_smooth(data = f,
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Ndep", y = "Ellenberg R") +
  scale_x_continuous(expand = c(0,0))
ggsave("Ndep effect on Ellenberg R score versions with data.png",
       path = "Outputs/Models/Difference/Univariate_measerror", 
       width = 16, height = 12, units = "cm")

# Read in models from files for pH plots
dir_mods <- "Outputs/Models/Difference/Univariate_measerror/"
mods <- list.files(path = dir_mods, pattern = "dep\\.rds$")
mods <- lapply(mods, function(x) readRDS(paste0(dir_mods, x)))

# remind myself which have divergent transitions
lapply(mods, function(x) {
  nuts_params(x) %>%
    filter(Parameter == "divergent__") %>%
    .$Value %>% sum()
})


pH_Ndep_mod <- mods[[1]]
pH_Sdep_mod <- mods[[2]]

summary(pH_Sdep_mod)
summary(pH_Ndep_mod)

mcmc_plot(pH_Sdep_mod, type = "acf")
mcmc_plot(pH_Ndep_mod, type = "acf")

sdep_pl <- plot(conditional_effects(pH_Sdep_mod), points = TRUE, plot = FALSE)[[3]]
ndep_pl <- plot(conditional_effects(pH_Ndep_mod), points = TRUE, plot = FALSE)[[3]]
sdep_pl + labs(x = expression("Change in S deposition (kg S ha"^-1*")"), 
               y = "pH change") + 
  scale_color_manual(values = mang_cols,
                     aesthetics = c("color","fill")) +
  ndep_pl + labs(x = expression("Cumulative N deposition (kg N ha"^-1*")"), 
                 y = "pH change")  + 
  scale_color_manual(values = mang_cols,
                     aesthetics = c("color","fill")) +
  plot_layout(guides = "collect")
ggsave("Sdep and Ndep effects on pH.png", path = dir_mods, 
       width = 15, height = 7, units = "cm", scale = 1.5)

plot_pars <- parnames(pH_Sdep_mod)[1:8]
param_summaries <- rbind(
  posterior_summary(pH_Sdep_mod) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    filter(Parameter %in% plot_pars) %>%
    mutate(Predictor = "Sdep"),
  posterior_summary(pH_Ndep_mod) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    filter(Parameter %in% plot_pars) %>%
    mutate(Predictor = "Ndep")
)
write.csv(param_summaries, paste0(dir_mods,"pH_NSdep_mods_params.csv"))

# Multivariate ####
# Combine pH and Ellenberg R with original pH and rain differences
ELL_pH <- PH %>%
  rename_with(~gsub("diff","PH_diff",.x)) %>%
  mutate(PHC_diff0719 = PHC_2019 - PHC_2007) %>%
  select(REP_ID, contains("diff7898"), contains("diff9807"),
         contains("diff0719")) %>%
  pivot_longer(contains("diff"), names_to = c("pH","Time_period"),
               names_sep= "_diff") %>%
  na.omit()  %>%
  pivot_wider(id_cols = c(REP_ID, Time_period), names_from = pH) %>%
  full_join(Ell_R) %>%
  mutate(Year = recode(Time_period,
                       "7898" = 1998,
                       "9807" = 2007,
                       "0719" = 2019),
         Year1 = recode(Time_period,
                        "7898" = 1978,
                        "9807" = 1998,
                        "0719" = 2007)) %>%
  left_join(BH_IMP) %>%
  left_join(filter(CS_plot_atdep, Duration == "8yr") %>%
              select(REP_ID, Year, Ndep, Sdep) %>%
              mutate(Year = ifelse(Year == 2018, 2019, Year))) %>%
  left_join(cs_rainfall_diff) %>%
  left_join(cs_rainfall_averages) %>%
  left_join(ELL_QA_diff) %>%
  left_join(PH_QA_diff) %>%
  left_join(rename(PH_long, Year1 = Year, Year1_pH = pH, 
                   Year1_pHCaCl2 = pH_CaCl2)) %>%
  left_join(rename(PH_long, Year2_pH = pH, 
                   Year2_pHCaCl2 = pH_CaCl2)) %>%
  left_join(CN)
summary(ELL_pH)
mice::md.pattern(ELL_pH)

# run model with rain differences
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W, ELL_SE = ELL_WH_W_SE) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, ELL_SE, Sdep, Ndep,
         fieldseason_rain, Year2_pH, PH, PH_DIW_SE, 
         N = NC_RATIO) %>%
  mutate(Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         N = as.numeric(scale(N)), 
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()
# 1237 obs

get_prior(bf(Ell | se(ELL_SE, sigma = TRUE) ~ me(PH,PH_DIW_SE,YRnm)*Year1_pH + (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID),
             family = "student") +
            bf(PH | se(PH_DIW_SE, sigma = TRUE)  ~ Sdep + fieldseason_rain + (1|p|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), family = "student") + 
            set_rescor(FALSE) + set_mecor(FALSE),
          data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(gamma(4,1), class = "nu", resp = "Ell"),
            prior(gamma(4,1), class = "nu", resp = "PH"),
            prior(normal(2,0.5), class = "ar", resp = "Ell"),
            prior(normal(2,0.5), class = "ar", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
            prior(lkj(2), class = "cor", group = "SQUARE"))

# prior predictive check
prior_mod <- brm(v, 
                 sample_prior = "only",
                 data = mod_data, prior = mod_pr, cores = 4, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50, resp = "Ell")
pp_check(prior_mod, nsamples = 50, resp = "PH")
plot(conditional_effects(prior_mod))
# underestimates peak at 0 on average but gets close
stancode(prior_mod) # check if what I think is happening is actually happening

# missing data model
make_stancode(bf(Ell | mi(ELL_SE) ~ mi(PH)*Year1_pH + (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID),
                 family = "student") +
                bf(PH | mi(PH_DIW_SE)  ~ Sdep + fieldseason_rain + (1|p|SQUARE) +
                     ar(time = YRnm, gr = REP_ID), family = "student") + 
                set_rescor(FALSE),
              data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(gamma(3,1), class = "nu", resp = "Ell"),
            prior(gamma(3,1), class = "nu", resp = "PH"),
            prior(normal(2,0.5), class = "ar", resp = "Ell"),
            prior(normal(2,0.5), class = "ar", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
            prior(lkj(2), class = "cor", group = "SQUARE"),
            prior(normal(0,1), class = "meanme", resp = "Ell"),
            prior(normal(0,1), class = "meanme", resp = "PH"),
            prior(student_t(3,0,3), class = "sdme", resp = "Ell"),
            prior(student_t(3,0,3), class = "sdme", resp = "PH"))

# prior predictive check
prior_mod <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH)*Year1_pH + (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID),
                    family = "student") +
                   bf(PH | mi(PH_DIW_SE)  ~ Sdep + fieldseason_rain + (1|p|SQUARE) +
                        ar(time = YRnm, gr = REP_ID), family = "student") + 
                   set_rescor(FALSE),
                 sample_prior = "only", save_all_pars = TRUE, save_mevars = TRUE,
                 data = mod_data, prior = mod_pr, cores = 4, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50, resp = "Ell")
pp_check(prior_mod, nsamples = 50, resp = "PH")
plot(conditional_effects(prior_mod))
# underestimates peak at 0 on average but gets close
stancode(prior_mod) # check if what I think is happening is actually happening

get_prior(bf(Ell | mi(ELL_SE) ~ mi(PH)*Year1_pH + (1|p|SQUARE),
             family = "student") +
            bf(PH | mi(PH_DIW_SE)  ~ Sdep + fieldseason_rain + (1|p|SQUARE), 
               family = "student") + 
            set_rescor(FALSE),
          data = mod_data)


# model run
phell_mod <- brm(bf(Ell | se(ELL_SE, sigma = TRUE) ~ PH*Year1_pH + (1|p|SQUARE) + 
                      ar(time = YRnm, gr = REP_ID, cov = TRUE), family = "student") +
                   bf(PH | se(PH_DIW_SE, sigma = TRUE) ~ Sdep + fieldseason_rain + (1|p|SQUARE) + 
                        ar(time = YRnm, gr = REP_ID, cov = TRUE), family = "student") +
                   set_rescor(FALSE),  save_all_pars = TRUE,
                 data = mod_data, prior = mod_pr, cores = 4, iter = 5000,
                 file = "Outputs/Models/Difference/WH_R_W/Ell_PH_multi_pHint_Sdepfieldrain")
summary(phell_mod)
plot(phell_mod)
pp_check(phell_mod, nsamples = 50, resp = "Ell")
pp_check(phell_mod, nsamples = 50, resp = "PH")

plot(conditional_effects(phell_mod), points = TRUE)

phell_mod <- add_criterion(phell_mod, "loo", reloo = TRUE, moment_match = TRUE,
                           file = "Outputs/Models/Difference/Ell_PH_multi_pHint_Sdepfieldrain")

# no year 1 pH effect
phell_mod_noint <- brm(bf(Ell ~ PH + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                          family = "student") +
                         bf(PH ~ Sdep + fieldseason_rain + (1|SQUARE) + 
                              ar(time = YRnm, gr = REP_ID), family = "student") +
                         set_rescor(FALSE),  save_all_pars = TRUE, 
                       data = mod_data, prior = mod_pr, cores = 4, iter = 5000,
                       file = "Outputs/Models/Difference/Ell_PH_multi_pH_Sdepfieldrain")
summary(phell_mod_noint)
plot(phell_mod_noint)
pp_check(phell_mod_noint, nsamples = 50, resp = "Ell")
pp_check(phell_mod_noint, nsamples = 50, resp = "PH")

plot(conditional_effects(phell_mod_noint), points = TRUE)

phell_mod_noint <- add_criterion(phell_mod_noint, "loo", reloo = TRUE, moment_match = TRUE,
                                 file = "Outputs/Models/Difference/Ell_PH_multi_pH_Sdepfieldrain")

loo_compare(phell_mod, phell_mod_noint)
# elpd_diff se_diff
# phell_mod        0.0       0.0   
# phell_mod_noint -0.6       2.2   

# ph model no rainfall
phell_mod_norain <- brm(bf(Ell ~ PH*Year1_pH + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                           family = "student") +
                          bf(PH ~ Sdep + (1|SQUARE) + 
                               ar(time = YRnm, gr = REP_ID), family = "student") +
                          set_rescor(FALSE),  save_all_pars = TRUE,
                        data = mod_data, prior = mod_pr, cores = 4, iter = 5000,
                        file = "Outputs/Models/Difference/Ell_PH_multi_pHint_Sdep")
summary(phell_mod_norain)
plot(phell_mod_norain)
pp_check(phell_mod_norain, nsamples = 50, resp = "Ell")
pp_check(phell_mod_norain, nsamples = 50, resp = "PH")

plot(conditional_effects(phell_mod_norain), points = TRUE)

phell_mod_norain <- add_criterion(phell_mod_norain, "loo", reloo = TRUE, moment_match = TRUE,
                                  file = "Outputs/Models/Difference/Ell_PH_multi_pHint_Sdep")


# ph model no rainfall no interaction
phell_mod_nointrain <- brm(bf(Ell ~ PH + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                              family = "student") +
                             bf(PH ~ Sdep + (1|SQUARE) + 
                                  ar(time = YRnm, gr = REP_ID), family = "student") +
                             set_rescor(FALSE),  save_all_pars = TRUE,
                           data = mod_data, prior = mod_pr, cores = 4, iter = 5000,
                           file = "Outputs/Models/Difference/Ell_PH_multi_pH_Sdep")
summary(phell_mod_nointrain)
plot(phell_mod_nointrain)
pp_check(phell_mod_nointrain, nsamples = 50, resp = "Ell")
pp_check(phell_mod_nointrain, nsamples = 50, resp = "PH")

plot(conditional_effects(phell_mod_nointrain), points = TRUE)

phell_mod_nointrain <- add_criterion(phell_mod_nointrain, "loo", reloo = TRUE, moment_match = TRUE,
                                     file = "Outputs/Models/Difference/Ell_PH_multi_pH_Sdep")

loo_compare(phell_mod, phell_mod_noint, phell_mod_nointrain, phell_mod_norain)
# elpd_diff se_diff
# phell_mod            0.0       0.0   
# phell_mod_noint     -0.6       2.2   
# phell_mod_norain    -3.4       3.1   
# phell_mod_nointrain -4.0       3.8  

# pH model no Sdep no interaction
phell_mod_nointsdep <- brm(bf(Ell ~ PH + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                              family = "student") +
                             bf(PH ~ fieldseason_rain + (1|SQUARE) + 
                                  ar(time = YRnm, gr = REP_ID), family = "student") +
                             set_rescor(FALSE),  save_all_pars = TRUE,
                           data = mod_data, prior = mod_pr, cores = 4, iter = 5000,
                           file = "Outputs/Models/Difference/Ell_PH_multi_pH_rain")
summary(phell_mod_nointsdep)
plot(phell_mod_nointsdep)
pp_check(phell_mod_nointsdep, nsamples = 50, resp = "Ell")
pp_check(phell_mod_nointsdep, nsamples = 50, resp = "PH")

plot(conditional_effects(phell_mod_nointsdep), points = TRUE)

phell_mod_nointsdep <- add_criterion(phell_mod_nointsdep, "loo", reloo = TRUE, moment_match = TRUE,
                                     file = "Outputs/Models/Difference/Ell_PH_multi_pH_rain")

loo_compare(phell_mod, phell_mod_noint, 
            phell_mod_nointrain, phell_mod_norain,
            phell_mod_nointsdep)
# elpd_diff se_diff
# phell_mod            0.0       0.0   
# phell_mod_nointsdep  0.0       2.3   
# phell_mod_noint     -0.6       2.2   
# phell_mod_norain    -3.4       3.1   
# phell_mod_nointrain -4.0       3.8 

loo_compare(phell_mod_nointrain, phell_mod_nointsdep)
# elpd_diff se_diff
# phell_mod_nointsdep  0.0       0.0   
# phell_mod_nointrain -4.0       3.0 

## adding Sdep and rainfall links to ellenberg R change
phell_mod_noint_sdep <- brm(bf(Ell ~ PH + Sdep + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                               family = "student") +
                              bf(PH ~ Sdep + fieldseason_rain + (1|SQUARE) + 
                                   ar(time = YRnm, gr = REP_ID), family = "student") +
                              set_rescor(FALSE),  save_all_pars = TRUE,
                            data = mod_data, prior = mod_pr, cores = 4, iter = 5000,
                            file = "Outputs/Models/Difference/Ell_PH_multi_EllSdep_pH_Sdepfieldrain")
summary(phell_mod_noint_sdep)
plot(phell_mod_noint_sdep)
pp_check(phell_mod_noint_sdep, nsamples = 50, resp = "Ell")
pp_check(phell_mod_noint_sdep, nsamples = 50, resp = "PH")

plot(conditional_effects(phell_mod_noint_sdep), points = TRUE)

phell_mod_noint_sdep <- add_criterion(phell_mod_noint_sdep, "loo", reloo = TRUE, moment_match = TRUE,
                                      file = "Outputs/Models/Difference/Ell_PH_multi_EllSdep_pH_Sdepfieldrain")

phell_mod_noint_rain <- brm(bf(Ell ~ PH + fieldseason_rain + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                               family = "student") +
                              bf(PH ~ Sdep + fieldseason_rain + (1|SQUARE) + 
                                   ar(time = YRnm, gr = REP_ID), family = "student") +
                              set_rescor(FALSE),  save_all_pars = TRUE,
                            data = mod_data, prior = mod_pr, cores = 4, iter = 5000,
                            file = "Outputs/Models/Difference/Ell_PH_multi_Ellrain_pH_Sdepfieldrain")
summary(phell_mod_noint_rain)
plot(phell_mod_noint_rain)
pp_check(phell_mod_noint_rain, nsamples = 50, resp = "Ell")
pp_check(phell_mod_noint_rain, nsamples = 50, resp = "PH")

plot(conditional_effects(phell_mod_noint_rain), points = TRUE)

phell_mod_noint_rain <- add_criterion(phell_mod_noint_rain, "loo", reloo = TRUE, moment_match = TRUE,
                                      file = "Outputs/Models/Difference/Ell_PH_multi_Ellrain_pH_Sdepfieldrain")


loo_compare(phell_mod_noint, phell_mod_noint_sdep, phell_mod_noint_rain)
#                       elpd_diff se_diff
# phell_mod_noint       0.0       0.0   
# phell_mod_noint_sdep -0.3       0.7   
# phell_mod_noint_rain -1.3       0.2 

# So Sdep and field season rainfall differences do not impact Ellenberg R
# directly (for whole plot weighted score)



# Multivariate with N #### 
# *** Low intensity management only ####
mod_data_low <- ELL_pH %>%
  filter(Management == "Low") %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH,
         ELL_SE = ELL_WH_W_SE_NORM, 
         PH_SE = PH_DIW_SE_NORM) %>%
  mutate(Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()

get_prior(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
               (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
            bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep, Ndep, Year1_pH) + 
                 (1|p|SQUARE) +
                 ar(time = YRnm, gr = REP_ID)) + 
            set_rescor(FALSE),
          data = mod_data_low)

mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(normal(2,0.5), class = "ar", resp = "Ell"),
            prior(normal(2,0.5), class = "ar", resp = "PH"),
            prior(student_t(3,0,1), class = "sd", resp = "Ell"),
            prior(student_t(3,0,1), class = "sd", resp = "PH"),
            prior(student_t(3,0,1), class = "sds", resp = "Ell"),
            prior(student_t(3,0,1), class = "sds", resp = "PH"),
            prior(student_t(3,0,1), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,1), class = "sigma", resp = "PH"),
            prior(student_t(3,0,5), class = "meanme", resp = "Ell"),
            prior(student_t(3,0,5), class = "meanme", resp = "PH"),
            prior(student_t(3,0,5), class = "sdme", resp = "Ell"),
            prior(student_t(3,0,5), class = "sdme", resp = "PH"),
            prior(lkj(2), class = "cor"))

full_mod <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                     (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                  bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep, Ndep, Year1_pH) + 
                       (1|p|SQUARE) +
                       ar(time = YRnm, gr = REP_ID)) + 
                  set_rescor(FALSE),
                data = mod_data_low, prior = mod_pr,
                cores = 4, iter = 4000, seed = 6856231456,
                file = "Outputs/Models/Difference/Test_fullmod_2903121")
# divergent transitions
summary(full_mod)
parnames(full_mod)
pairs(full_mod$fit, pars = c("b_Ell_Intercept",
                             "b_PH_Intercept",
                             "sd_SQUARE__PH_Intercept",
                             "sd_SQUARE__Ell_Intercept",
                             "ar_PH[1]",
                             "sds_PH_sSdepNdepYear1_pH_1",
                             "sigma_PH",
                             "bs_PH_sSdepNdepYear1_pH_9",
                             "b_PH_fieldseason_rain"),
      log = TRUE)
plot(full_mod, pars = grep("^s_PH", parnames(full_mod), value = TRUE),
     fixed = TRUE)
plot(full_mod2, pars = parnames(full_mod2),
     window = c(0,750), chain = 3)

partition_fullmod <- util$partition_div(full_mod$fit)

spline_vars <- grep("^s_PH_sSdepNdepYear1_pH_1",
                    colnames(partition_fullmod[[1]]),
                    value = TRUE)
par(mfrow=c(5,5))
for (i in 1:100){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1], 
       log(partition_fullmod[[2]]$sds_PH_sSdepNdepYear1_pH_1),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_pH_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_PH_sSdepNdepYear1_pH_1),
         col=c_yellow_trans, pch=16)
}

spline_vars <- grep("^s_Ell_sNdepYear1_pH_",
                    parnames(full_mod3),
                    value = TRUE)
par(mfrow=c(3,3), mar = c(4,4,2,2))
for (i in 1:27){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1],
       log(partition_fullmod[[2]][,"sds_Ell_sNdepYear1_pH_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_Ell_sNdepYear1_pH_1),
         col=c_yellow_trans, pch=16)
}

# plot MANY MUCH THING
plot_pars <- parnames(full_mod)[1:24]

pairs(full_mod$fit, pars = plot_pars[1:12],
      log = TRUE)
pairs(full_mod$fit, pars = plot_pars[13:24],
      log = TRUE)


par(mfrow=c(3,4))
for(i in 1:12){
  name_x <- plot_pars[i]
  for(j in 13:24){
    name_y <- plot_pars[j]
    y <- if(grepl("^s",name_y))
      log(partition_fullmod[[2]][name_y][,1]) else
                partition_fullmod[[2]][name_y][,1]
    ydiv <- if(grepl("^s",name_y)) 
      log(partition_fullmod[[1]][name_y][,1]) else
        partition_fullmod[[1]][name_y][,1]
    
    plot(partition_fullmod[[2]][name_x][,1],
         y, xlab = name_x, ylab = name_y,
         col=c_dark_trans, pch=16, main="")
    points(partition_fullmod[[1]][name_x][,1],
           ydiv,
           col=c_yellow_trans, pch=16)
  }
  invisible(readline(prompt="Press [enter] to continue"))
}

plot_pars <- parnames(full_mod)
par(mfrow=c(5,5), mar = c(4,4,1,1))
for(i in 4:24){
  name_x <- plot_pars[i]
  x <- if(grepl("^s",name_x))
    log(partition_fullmod[[2]][name_x][,1]) else
      partition_fullmod[[2]][name_x][,1]
  xdiv <- if(grepl("^s",name_x)) 
    log(partition_fullmod[[1]][name_x][,1]) else
      partition_fullmod[[1]][name_x][,1]
  
  for(j in 350:499){
    name_y <- plot_pars[j]
    y <- partition_fullmod[[2]][name_y][,1]
    ydiv <- partition_fullmod[[1]][name_y][,1]
    
    plot(x,
         y, xlab = name_x, ylab = name_y,
         col=c_dark_trans, pch=16, main="")
    points(xdiv,
           ydiv,
           col=c_yellow_trans, pch=16)
    
    if(j %in% c(5*5*(1:19)+24))
      invisible(readline(prompt="Press [enter] to continue"))
  }
  
}


full_mod2 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                      (1|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                   bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep, Ndep, Year1_pH) + 
                        (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID)) + 
                   set_rescor(FALSE),
                 data = mod_data_low, prior = mod_pr[1:16,],
                 cores = 4, iter = 4000, seed = 4575421286,
                 file = "Outputs/Models/Difference/Test_fullmod2_2903121")
# divergent transitions
summary(full_mod2)
parnames(full_mod2)
pairs(full_mod2$fit, pars = c("b_Ell_Intercept",
                              "b_PH_Intercept",
                              "sd_SQUARE__PH_Intercept",
                              "sd_SQUARE__Ell_Intercept",
                              "ar_PH[1]",
                              "sds_PH_sSdepNdepYear1_pH_1",
                              "sigma_PH",
                              "bs_PH_sSdepNdepYear1_pH_9",
                              "b_PH_fieldseason_rain"),
      log = TRUE)
plot(full_mod2, pars = c("b_PH_Intercept",
                         "sd_SQUARE__PH_Intercept",
                         "ar_PH[1]",
                         "sds_PH_sSdepNdepYear1_pH_1",
                         "sigma_PH",
                         "bs_PH_sSdepNdepYear1_pH_9",
                         "b_PH_fieldseason_rain"))
plot(full_mod2, pars = parnames(full_mod2),
     window = c(0,750), chain = 3)


partition_fullmod <- util$partition_div(full_mod2$fit)

spline_vars <- grep("^s_PH_sSdepNdepYear1_pH_1",
                    colnames(partition_fullmod[[1]]),
                    value = TRUE)
par(mfrow=c(5,5), mar = c(4,4,2,2))
for (i in 1:100){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1], 
       log(partition_fullmod[[2]]$sds_PH_sSdepNdepYear1_pH_1),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_pH_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_PH_sSdepNdepYear1_pH_1),
         col=c_green_trans, pch=16)
}


mod2 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                 (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
              bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep) + s(Ndep, Year1_pH) + 
                   (1|p|SQUARE) +
                   ar(time = YRnm, gr = REP_ID)) + 
              set_rescor(FALSE),
            data = mod_data_low, prior = mod_pr,
            cores = 4, iter = 4000, seed = 6856231456,
            file = "Outputs/Models/Difference/Test_mod2_2903121")
# divergent transitions
summary(mod2)
parnames(mod2)
pairs(mod2$fit, pars = c("b_Ell_Intercept",
                         "b_PH_Intercept",
                         "sd_SQUARE__PH_Intercept",
                         "sd_SQUARE__Ell_Intercept",
                         "ar_PH[1]",
                         "sds_PH_sNdepYear1_pH_1",
                         "sigma_PH",
                         "bs_PH_sNdepYear1_pH_2",
                         "b_PH_fieldseason_rain"),
      log = TRUE)
plot(mod2, pars = c("b_PH_Intercept",
                    "sd_SQUARE__PH_Intercept",
                    "ar_PH[1]",
                    "sds_PH_sNdepYear1_pH_1",
                    "sigma_PH",
                    "bs_PH_sNdepYear1_pH_2",
                    "b_PH_fieldseason_rain"))
plot(mod2, pars = parnames(mod2), fixed = TRUE,
     window = c(500,1500))
pairs(mod2$fit, pars = c("sds_PH_sSdep_1",
                         "bs_PH_sSdep_1",
                         "ar_PH[1]",
                         "b_PH_fieldseason_rain",
                         "sd_SQUARE__PH_Intercept",
                         "r_SQUARE__PH[63,Intercept]",
                         "s_PH_sNdepYear1_pH_1[8]",
                         "s_PH_sSdep_1[7]"),
      log = TRUE)
pairs(mod2$fit, pars = c("sds_PH_sSdep_1",
                         "s_PH_sSdep_1[1]",
                         "s_PH_sSdep_1[2]",
                         "s_PH_sSdep_1[3]",
                         "s_PH_sSdep_1[4]",
                         "s_PH_sSdep_1[5]",
                         "s_PH_sSdep_1[6]",
                         "s_PH_sSdep_1[7]",
                         "s_PH_sSdep_1[8]"),
      log = TRUE)

mod_data_low %>%
  mutate(Weird = ifelse(SQUARE %in% c("1057","1152","161","215","203",
                                      "205","459","467","63","820","877"),
                        1,0)) %>%
  ggplot(aes(x = Sdep, y = PH, colour = Weird)) + 
  geom_point()


partition_mod2 <- util$partition_div(mod2$fit)

spline_vars <- grep("^s_PH_sSdep_1",
                    colnames(partition_mod2[[1]]),
                    value = TRUE)
par(mfrow=c(2,4))
for (i in 1:length(spline_vars)){
  name_x <- spline_vars[i]
  plot(partition_mod2[[2]][name_x][,1], 
       log(partition_mod2[[2]]$sds_PH_sSdep_1),
       col=c_dark_trans, pch=16, main="",
       xlim = c(-7,7), ylim = c(-10,2),
       ylab="log(sds_Sdep)", xlab = paste("Sdep spline:",i))
  points(partition_mod2[[1]][name_x][,1], 
         log(partition_mod2[[1]]$sds_PH_sSdep_1),
         col=c_yellow_trans, pch=16)
}

spline_vars <- grep("^s_PH_sNdepYear1_pH_1",
                    colnames(partition_mod2[[1]]),
                    value = TRUE)
par(mfrow=c(3,3))
for (i in 1:length(spline_vars)){
  name_x <- spline_vars[i]
  plot(partition_mod2[[2]][name_x][,1], 
       log(partition_mod2[[2]]$sds_PH_sSdep_1),
       col=c_dark_trans, pch=16, main="",
       # xlim = c(-7,7), ylim = c(-10,2),
       ylab="log(sds_NdeppH)", xlab = paste("Ndep pH spline:",i))
  points(partition_mod2[[1]][name_x][,1], 
         log(partition_mod2[[1]]$sds_PH_sSdep_1),
         col=c_yellow_trans, pch=16)
}

spline_vars <- grep("^s_Ell_sNdepYear1_pH_",
                    parnames(mod2),
                    value = TRUE)
par(mfrow=c(3,3), mar = c(4,4,2,2))
for (i in 1:27){
  name_x <- spline_vars[i]
  plot(partition_mod2[[2]][name_x][,1],
       log(partition_mod2[[2]][,"sds_Ell_sNdepYear1_pH_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
  points(partition_mod2[[1]][name_x][,1], 
         log(partition_mod2[[1]]$sds_Ell_sNdepYear1_pH_1),
         col=c_yellow_trans, pch=16)
}




mod3 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                 (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
              bf(PH | mi(PH_SE)  ~ fieldseason_rain + Sdep + s(Ndep, Year1_pH) + 
                   (1|p|SQUARE) +
                   ar(time = YRnm, gr = REP_ID)) + 
              set_rescor(FALSE),
            data = mod_data_low, prior = mod_pr,
            cores = 4, iter = 4000, seed = 6856231456,
            file = "Outputs/Models/Difference/Test_mod3_2903121")
summary(mod3)
plot(mod3)

mod3_pr <- c(mod_pr[1:16,], 
             prior(lkj(5), class = "cor"))
mod3 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                 (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
              bf(PH | mi(PH_SE)  ~ fieldseason_rain + Sdep + s(Ndep, Year1_pH) + 
                   (1|p|SQUARE) +
                   ar(time = YRnm, gr = REP_ID)) + 
              set_rescor(FALSE),
            data = mod_data_low, prior = mod3_pr,
            cores = 4, iter = 4000, seed = 6856231456,
            file = "Outputs/Models/Difference/Test_mod3pr_2903121")
summary(mod3)
plot(mod3)

mod3v2 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                   (1|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                bf(PH | mi(PH_SE)  ~ fieldseason_rain + Sdep + s(Ndep, Year1_pH) + 
                     (1|SQUARE) +
                     ar(time = YRnm, gr = REP_ID)) + 
                set_rescor(FALSE),
              data = mod_data_low, prior = mod_pr[1:16,],
              cores = 4, iter = 4000, seed = 6856231456,
              file = "Outputs/Models/Difference/Test_mod3v2_2903121")
summary(mod3v2)
plot(mod3v2)
plot(conditional_effects(mod3v2))

spline_vars <- grep("^s_Ell_sNdepYear1_pH_",
                    parnames(mod3v2),
                    value = TRUE)
partition <- util$partition_div(mod3v2$fit)
par(mfrow=c(3,3), mar = c(4,4,2,2))
for (i in 1:27){
  name_x <- spline_vars[i]
  plot(partition[[2]][name_x][,1],
       log(partition[[2]][,"sds_Ell_sNdepYear1_pH_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
  points(partition[[1]][name_x][,1], 
         log(partition[[1]]$sds_Ell_sNdepYear1_pH_1),
         col=c_yellow_trans, pch=16)
}



mod_pr2 <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
             prior(normal(0,0.5), class = "b", resp = "PH"),
             prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
             prior(normal(0,0.25), class = "Intercept", resp = "PH"),
             prior(normal(2,0.5), class = "ar", resp = "Ell"),
             prior(normal(2,0.5), class = "ar", resp = "PH"),
             prior(student_t(3,0,1), class = "sd", resp = "Ell"),
             prior(student_t(3,0,1), class = "sd", resp = "PH"),
             prior(normal(0,1), class = "sds", resp = "Ell"),
             prior(normal(0,1), class = "sds", resp = "PH"),
             prior(student_t(3,0,1), class = "sigma", resp = "Ell"),
             prior(student_t(3,0,1), class = "sigma", resp = "PH"),
             prior(student_t(3,0,5), class = "meanme", resp = "Ell"),
             prior(student_t(3,0,5), class = "meanme", resp = "PH"),
             prior(student_t(3,0,5), class = "sdme", resp = "Ell"),
             prior(student_t(3,0,5), class = "sdme", resp = "PH"),
             prior(lkj(5), class = "cor"))
mod2v2 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                   (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep) + s(Ndep, Year1_pH) + 
                     (1|p|SQUARE) +
                     ar(time = YRnm, gr = REP_ID)) + 
                set_rescor(FALSE),
              data = mod_data_low, prior = mod_pr2,
              cores = 4, iter = 4000, seed = 6856231456,
              file = "Outputs/Models/Difference/Test_mod2v2_2903121")
# divergent transitions
summary(mod2v2)
parnames(mod2v2)
pairs(mod2v2$fit, pars = c("b_Ell_Intercept",
                           "b_PH_Intercept",
                           "sd_SQUARE__PH_Intercept",
                           "sd_SQUARE__Ell_Intercept",
                           "ar_PH[1]",
                           "sds_PH_sNdepYear1_pH_1",
                           "sigma_PH",
                           "bs_PH_sNdepYear1_pH_2",
                           "b_PH_fieldseason_rain"),
      log = TRUE)
plot(mod2, pars = c("b_PH_Intercept",
                    "sd_SQUARE__PH_Intercept",
                    "ar_PH[1]",
                    "sds_PH_sNdepYear1_pH_1",
                    "sigma_PH",
                    "bs_PH_sNdepYear1_pH_2",
                    "b_PH_fieldseason_rain"))
plot(mod2v2, pars = parnames(mod2), fixed = TRUE,
     window = c(500,1500))
pairs(mod2v2$fit, pars = c("sds_PH_sSdep_1",
                           "bs_PH_sSdep_1",
                           "ar_PH[1]",
                           "b_PH_fieldseason_rain",
                           "sd_SQUARE__PH_Intercept",
                           "r_SQUARE__PH[63,Intercept]",
                           "s_PH_sNdepYear1_pH_1[8]",
                           "s_PH_sSdep_1[7]"),
      log = TRUE)
pairs(mod2v2$fit, pars = c("sds_PH_sSdep_1",
                           "s_PH_sSdep_1[1]",
                           "s_PH_sSdep_1[2]",
                           "s_PH_sSdep_1[3]",
                           "s_PH_sSdep_1[4]",
                           "s_PH_sSdep_1[5]",
                           "s_PH_sSdep_1[6]",
                           "s_PH_sSdep_1[7]",
                           "s_PH_sSdep_1[8]"),
      log = TRUE)

mod_data_low %>%
  mutate(Weird = ifelse(SQUARE %in% c("1057","1152","161","215","203",
                                      "205","459","467","63","820","877"),
                        1,0)) %>%
  ggplot(aes(x = Sdep, y = PH, colour = Weird)) + 
  geom_point()

mod4 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                 (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
              bf(PH | mi(PH_SE)  ~ s(Sdep) + 
                   (1|p|SQUARE) +
                   ar(time = YRnm, gr = REP_ID)) + 
              set_rescor(FALSE),
            data = mod_data_low, prior = mod3_pr,
            cores = 4, iter = 4000, seed = 6856231456,
            file = "Outputs/Models/Difference/Test_mod4_2903121")


partition_mod <- util$partition_div(mod4$fit)

spline_vars <- grep("^s_PH_sSdep_1",
                    colnames(partition_mod[[1]]),
                    value = TRUE)
par(mfrow=c(2,4))
for (i in 1:length(spline_vars)){
  name_x <- spline_vars[i]
  plot(partition_mod[[2]][name_x][,1], 
       log(partition_mod[[2]]$sds_PH_sSdep_1),
       col=c_dark_trans, pch=16, main="",
       xlim = c(-7,7), ylim = c(-10,2),
       ylab="log(sds_Sdep)", xlab = paste("Sdep spline:",i))
  points(partition_mod[[1]][name_x][,1], 
         log(partition_mod[[1]]$sds_PH_sSdep_1),
         col=c_yellow_trans, pch=16)
}

mod5 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                 (1|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
              bf(PH | mi(PH_SE)  ~ s(Sdep) ) + 
              set_rescor(FALSE),
            data = mod_data_low, prior = mod3_pr[c(1:5,7,9,11:16),],
            cores = 4, iter = 4000, seed = 6856231456)

mod6 <- brm(PH | mi(PH_SE)  ~ s(Sdep),
            data = mod_data_low, 
            prior = prior(normal(0,1), class = "b"),
            cores = 4, iter = 4000, seed = 6856231456)
pairs(mod6$fit, log = TRUE)

mod7 <- brm(PH ~ s(Sdep),
            data = mod_data_low, 
            prior = prior(normal(0,1), class = "b"),
            cores = 4, iter = 4000, seed = 6856231456)

full_mod3 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Ndep, Year1_pH) + 
                      (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                   bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep, Ndep, Year1_pH) + 
                        (1|p|SQUARE) +
                        ar(time = YRnm, gr = REP_ID)) + 
                   set_rescor(FALSE),
                 data = mod_data_low, prior = mod3_pr,
                 control = list(adapt_delta = 0.99999999, max_treedepth = 20),
                 cores = 4, iter = 4000, seed = 6856231456,
                 file = "Outputs/Models/Difference/Test_fullmod3_3103121")

spline_vars <- grep("^s_PH_sSdepNdepYear1_pH_1",
                    parnames(full_mod3),
                    value = TRUE)
test <- as.matrix(full_mod3)
par(mfrow=c(4,5), mar = c(4,4,2,2))
for (i in 1:100){
  name_x <- spline_vars[i]
  plot(test[,name_x], 
       log(test[,"sds_PH_sSdepNdepYear1_pH_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_pH_spline)", xlab = paste("Spline:",i))
}

spline_vars <- grep("^bs_PH_sSdepNdepYear1_pH_",
                    parnames(full_mod3),
                    value = TRUE)
par(mfrow=c(3,3), mar = c(4,4,2,2))
for (i in 1:9){
  name_x <- spline_vars[i]
  plot(test[,name_x], 
       log(test[,"sds_PH_sSdepNdepYear1_pH_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_pH_spline)", xlab = paste("Spline:",i))
}

spline_vars <- grep("^s_Ell_sNdepYear1_pH_",
                    parnames(full_mod3),
                    value = TRUE)
par(mfrow=c(3,3), mar = c(4,4,2,2))
for (i in 1:27){
  name_x <- spline_vars[i]
  plot(test[,name_x], 
       log(test[,"sds_Ell_sNdepYear1_pH_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
}


# 3d surface
library(plotly)

nd <-
  tibble(Year1_pH = rep(c(-1,0,1), each = 400),
         Ndep = seq(-1.7,3.13, length.out = 20) %>% 
           rep(each = 20) %>%
           rep(times = 3),
         Sdep = seq(-6.13,1.25, length.out = 20) %>% 
           rep(times = 20) %>%
           rep(times = 3),
         REP_ID = 1:1200,
         fieldseason_rain = 0,
         ELL_SE = 0.3, PH_SE = 0.3,
         PH = 0, YRnm = 2)
mod_fit <-  fitted(full_mod3, newdata = nd, re_formula = NA) %>%
  as_tibble() %>%
  bind_cols(nd) 

mod_fit <- list(Ndep = seq(-1.7,3.13, length.out = 20),
                Sdep = seq(-6.13,1.25, length.out = 20),
                Est_lowpH = matrix(mod_fit$Estimate.PH[1:400],
                                   nrow = 20, ncol= 20),
                Est_medpH = matrix(mod_fit$Estimate.PH[401:800],
                                   nrow = 20, ncol= 20),
                Est_highpH = matrix(mod_fit$Estimate.PH[801:1200],
                                    nrow = 20, ncol= 20),
                Q025_lowpH = matrix(mod_fit$Q2.5.PH[1:400],
                                    nrow = 20, ncol= 20),
                Q025_medpH = matrix(mod_fit$Q2.5.PH[401:800],
                                    nrow = 20, ncol= 20),
                Q025_highpH = matrix(mod_fit$Q2.5.PH[801:1200],
                                     nrow = 20, ncol= 20),
                Q975_lowpH = matrix(mod_fit$Q97.5.PH[1:400],
                                    nrow = 20, ncol= 20),
                Q975_medpH = matrix(mod_fit$Q97.5.PH[401:800],
                                    nrow = 20, ncol= 20),
                Q975_highpH = matrix(mod_fit$Q97.5.PH[801:1200],
                                     nrow = 20, ncol= 20))
# no CI
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(x = mod_fit$Ndep,
                           y = mod_fit$Sdep,
                           z = mod_fit$Est_lowpH,
                           name = "Low pH",
                           colorscale = list(c(0,1),c("rgb(240, 123, 31)","rgb(196, 86, 0)")))

fig <- fig %>% add_surface(x = mod_fit$Ndep,
                           y = mod_fit$Sdep,
                           z = mod_fit$Est_medpH,
                           name = "Mean pH",
                           colorscale = list(c(0,1),c("rgb(9, 237, 175","rgb(0, 156, 113)")))

fig <- fig %>% add_surface(x = mod_fit$Ndep,
                           y = mod_fit$Sdep,
                           z = mod_fit$Est_highpH,
                           name = "High pH",
                           colorscale = list(c(0,1), c("rgb(0, 145, 227)","rgb(2, 98, 153)")))

fig <- fig %>% add_markers(x = filter(mod_data_low, Year1_pH < -0.5)$Ndep,
                           y = filter(mod_data_low, Year1_pH < -0.5)$Sdep,
                           z = filter(mod_data_low, Year1_pH < -0.5)$PH,
                           name = "Low pH",
                           color = "#D55E00",
                           showlegend = FALSE)

fig <- fig %>% add_markers(x = filter(mod_data_low, Year1_pH < 0.5 & Year1_pH > -0.5)$Ndep,
                           y = filter(mod_data_low, Year1_pH < 0.5 & Year1_pH > -0.5)$Sdep,
                           z = filter(mod_data_low, Year1_pH < 0.5 & Year1_pH > -0.5)$PH,
                           name = "Mean pH",
                           color = "#009E73",
                           showlegend = FALSE)

fig <- fig %>% add_markers(x = filter(mod_data_low, Year1_pH > 0.5)$Ndep,
                           y = filter(mod_data_low, Year1_pH > 0.5)$Sdep,
                           z = filter(mod_data_low, Year1_pH > 0.5)$PH,
                           name = "High pH",
                           color = "#0072B2",
                           showlegend = FALSE)

fig <- fig %>% layout(scene = list(xaxis = list(title = "Ndep"), 
                                   yaxis = list(title = "Sdep"),
                                   zaxis = list(title = "pH change")))
fig

# with CI
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_markers(x = mod_data_low$Ndep,
                           y = mod_data_low$Sdep,
                           z = mod_data_low$PH,
                           name = "Low pH",
                           color = mod_data_low$Year1_pH<0,
                           colors = c("#0072B2","#D55E00"),
                           opacity = 0.5,
                           showlegend = FALSE)

fig <- fig %>% add_surface(x = mod_fit$Ndep,
                           y = mod_fit$Sdep,
                           z = mod_fit$Est_lowpH,
                           colorscale = list(c(0,1),c("rgb(240, 123, 31)","rgb(196, 86, 0)")))
# colors = "#D55E00")
fig <- fig %>% add_surface(x = mod_fit$Ndep,
                           y = mod_fit$Sdep,
                           z = mod_fit$Q025_lowpH, opacity = 0.25,
                           colorscale = list(c(0,1),c("rgb(240, 123, 31)","rgb(196, 86, 0)")))
fig <- fig %>% add_surface(x = mod_fit$Ndep,
                           y = mod_fit$Sdep,
                           z = mod_fit$Q975_lowpH, opacity = 0.25,
                           colorscale = list(c(0,1),c("rgb(240, 123, 31)","rgb(196, 86, 0)")))

fig <- fig %>% add_surface(x = mod_fit$Ndep,
                           y = mod_fit$Sdep,
                           z = mod_fit$Est_highpH,
                           colorscale = list(c(0,1), c("rgb(0, 145, 227)","rgb(2, 98, 153)")))
fig <- fig %>% add_surface(x = mod_fit$Ndep,
                           y = mod_fit$Sdep,
                           z = mod_fit$Q025_highpH, opacity = 0.25,
                           colorscale = list(c(0,1), c("rgb(0, 145, 227)","rgb(2, 98, 153)")))
fig <- fig %>% add_surface(x = mod_fit$Ndep,
                           y = mod_fit$Sdep,
                           z = mod_fit$Q975_highpH, opacity = 0.25,
                           colorscale = list(c(0,1), c("rgb(0, 145, 227)","rgb(2, 98, 153)")))


fig <- fig %>% layout(scene = list(xaxis = list(title = "Ndep"), 
                                   yaxis = list(title = "Sdep"),
                                   zaxis = list(title = "pH change")))
fig


mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         ELL_SE = ELL_WH_W_SE_NORM, 
         PH_SE = PH_DIW_SE_NORM) %>%
  mutate(Management = as.factor(Management),
         Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()

mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(normal(2,0.5), class = "ar", resp = "Ell"),
            prior(normal(2,0.5), class = "ar", resp = "PH"),
            prior(student_t(3,0,1), class = "sd", resp = "Ell"),
            prior(student_t(3,0,1), class = "sd", resp = "PH"),
            prior(student_t(3,0,1), class = "sds", resp = "Ell"),
            prior(student_t(3,0,1), class = "sds", resp = "PH"),
            prior(student_t(3,0,1), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,1), class = "sigma", resp = "PH"),
            prior(student_t(3,0,5), class = "meanme", resp = "Ell"),
            prior(student_t(3,0,5), class = "meanme", resp = "PH"),
            prior(student_t(3,0,5), class = "sdme", resp = "Ell"),
            prior(student_t(3,0,5), class = "sdme", resp = "PH"),
            prior(lkj(5), class = "cor"))

full_mod_hab <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH)*Management + s(Ndep, Year1_pH, by = Management) + 
                     (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                  bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep, Ndep, Year1_pH, by = Management) + 
                       (1|p|SQUARE) +
                       ar(time = YRnm, gr = REP_ID)) + 
                  set_rescor(FALSE),
                data = mod_data, prior = mod_pr,
                cores = 4, iter = 4000, seed = 6856231456,
                file = "Outputs/Models/Difference/Test_fullmodhabv2_010421")
summary(full_mod_hab)

partition_fullmod <- util$partition_div(full_mod_hab$fit)

spline_vars <- grep("^s_PH_sSdepNdepYear1_pH",
                    colnames(partition_fullmod[[1]]),
                    value = TRUE)
par(mfrow=c(5,5))
for (i in 1:100){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1], 
       log(partition_fullmod[[2]]$sds_PH_sSdepNdepYear1_pHManagementHigh_1),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_pH_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_PH_sSdepNdepYear1_pHManagementHigh_1),
         col=c_yellow_trans, pch=16)
}
for (i in 101:200){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1], 
       log(partition_fullmod[[2]]$sds_PH_sSdepNdepYear1_pHManagementLow_1),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_pH_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_PH_sSdepNdepYear1_pHManagementLow_1),
         col=c_yellow_trans, pch=16)
}

spline_vars <- grep("^s_Ell_sNdepYear1_pH",
                    parnames(full_mod_hab),
                    value = TRUE)
par(mfrow=c(3,3))
for (i in 1:27){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1],
       log(partition_fullmod[[2]][,"sds_Ell_sNdepYear1_pHManagementHigh_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_Ell_sNdepYear1_pHManagementHigh_1),
         col=c_yellow_trans, pch=16)
}
for (i in 28:54){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1],
       log(partition_fullmod[[2]][,"sds_Ell_sNdepYear1_pHManagementLow_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_Ell_sNdepYear1_pHManagementLow_1),
         col=c_yellow_trans, pch=16)
}


mod3_hab <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH)*Management + s(Ndep, Year1_pH, by = Management) + 
                         (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                      bf(PH | mi(PH_SE)  ~ fieldseason_rain + Sdep*Management + s(Ndep, Year1_pH, by = Management) + 
                           (1|p|SQUARE) +
                           ar(time = YRnm, gr = REP_ID)) + 
                      set_rescor(FALSE),
                    data = mod_data, prior = mod_pr,
                    cores = 4, iter = 4000, seed = 6856231456,
                    file = "Outputs/Models/Difference/Test_mod3hab_010421")
summary(mod3_hab)
mod3_stepsizes <- sapply(1:4, function(c) 
  get_sampler_params(mod3_hab$fit, inc_warmup=FALSE)[[c]][,'stepsize__'][1])
partition_fullmod <- util$partition_div(mod3_hab$fit)

spline_vars <- grep("^s_Ell_sNdepYear1_pH",
                    parnames(mod3_hab),
                    value = TRUE)
par(mfrow=c(3,3))
for (i in 1:27){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1],
       log(partition_fullmod[[2]][,"sds_Ell_sNdepYear1_pHManagementHigh_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_Ell_sNdepYear1_pHManagementHigh_1),
         col=c_yellow_trans, pch=16)
}
for (i in 28:54){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1],
       log(partition_fullmod[[2]][,"sds_Ell_sNdepYear1_pHManagementLow_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_Ell_sNdepYear1_pHManagementLow_1),
         col=c_yellow_trans, pch=16)
}



mod_pr2 <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(normal(2,0.5), class = "ar", resp = "Ell"),
            prior(normal(2,0.5), class = "ar", resp = "PH"),
            prior(student_t(3,0,1), class = "sd", resp = "Ell"),
            prior(student_t(3,0,1), class = "sd", resp = "PH"),
            prior(gamma(3,2), class = "sds", resp = "Ell"),
            prior(student_t(3,0,1), class = "sds", resp = "PH"),
            prior(student_t(3,0,1), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,1), class = "sigma", resp = "PH"),
            prior(student_t(3,0,5), class = "meanme", resp = "Ell"),
            prior(student_t(3,0,5), class = "meanme", resp = "PH"),
            prior(student_t(3,0,5), class = "sdme", resp = "Ell"),
            prior(student_t(3,0,5), class = "sdme", resp = "PH"),
            prior(lkj(5), class = "cor"))

full_mod_hab_v2 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH)*Management + s(Ndep, Year1_pH, by = Management) + 
                            (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                         bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep, Ndep, Year1_pH, by = Management) + 
                              (1|p|SQUARE) +
                              ar(time = YRnm, gr = REP_ID)) + 
                         set_rescor(FALSE),
                       data = mod_data, prior = mod_pr,
                       cores = 4, iter = 4000, seed = 6856231456,
                       file = "Outputs/Models/Difference/Test_fullmodhabv2_010421")
summary(full_mod_hab_v2)
plot(conditional_smooths(full_mod_hab_v2, too_far = 0.1))

partition_fullmod <- util$partition_div(full_mod_hab_v2$fit)

spline_vars <- grep("^s_Ell_sNdepYear1_pH",
                    parnames(full_mod_hab_v2),
                    value = TRUE)
par(mfrow=c(3,3))
for (i in 1:27){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1],
       log(partition_fullmod[[2]][,"sds_Ell_sNdepYear1_pHManagementHigh_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_Ell_sNdepYear1_pHManagementHigh_1),
         col=c_yellow_trans, pch=16)
}
for (i in 28:54){
  name_x <- spline_vars[i]
  plot(partition_fullmod[[2]][name_x][,1],
       log(partition_fullmod[[2]][,"sds_Ell_sNdepYear1_pHManagementLow_1"]),
       col=c_dark_trans, pch=16, main="",
       ylab="log(sds_Ell_spline)", xlab = paste("Spline:",i))
  points(partition_fullmod[[1]][name_x][,1], 
         log(partition_fullmod[[1]]$sds_Ell_sNdepYear1_pHManagementLow_1),
         col=c_yellow_trans, pch=16)
}

full_mod_hab_v3 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH)*Management + 
                            (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                         bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep, Ndep, Year1_pH, by = Management) + 
                              (1|p|SQUARE) +
                              ar(time = YRnm, gr = REP_ID)) + 
                         set_rescor(FALSE),
                       data = mod_data, prior = mod_pr2[-9,],
                       cores = 4, iter = 5000, seed = 6856231456,
                       file = "Outputs/Models/Difference/Test_fullmodhabv3_130421")
summary(full_mod_hab_v3)
partition_fullmod <- util$partition_div(full_mod_hab_v3$fit)

spline_vars <- grep("^s_PH_s",
                    parnames(full_mod_hab_v3),
                    value = TRUE)
par(mfrow=c(4,5))
for(i in 1:200){
  name_x <- spline_vars[i]
  name_y <- ifelse(i <=100,"sds_PH_sSdepNdepYear1_pHManagementHigh_1",
                   "sds_PH_sSdepNdepYear1_pHManagementLow_1")
  y <- log(partition_fullmod[[2]][name_y][,1]) 
  ydiv <- log(partition_fullmod[[1]][name_y][,1]) 
  
  plot(partition_fullmod[[2]][name_x][,1],
       y, xlab = name_x, ylab = name_y,
       col=c_dark_trans, pch=16, main="")
  points(partition_fullmod[[1]][name_x][,1],
         ydiv,
         col=c_yellow_trans, pch=16)
  if(i %in% c(4*5*(1:9)))
  invisible(readline(prompt="Press [enter] to continue"))
}

full_mod_hab_v4 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH)*Management + 
                            (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                         bf(PH | mi(PH_SE)  ~ fieldseason_rain + s(Sdep, Ndep, by = Management) + 
                              (1|p|SQUARE) +
                              ar(time = YRnm, gr = REP_ID)) + 
                         set_rescor(FALSE),
                       data = mod_data, prior = mod_pr2[-9,],
                       cores = 4, iter = 5000, seed = 6856231456,
                       file = "Outputs/Models/Difference/Test_fullmodhabv4_130421")
summary(full_mod_hab_v4)


# linear models
mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
             prior(normal(0,0.5), class = "b", resp = "PH"),
             prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
             prior(normal(0,0.25), class = "Intercept", resp = "PH"),
             prior(normal(2,0.5), class = "ar", resp = "Ell"),
             prior(normal(2,0.5), class = "ar", resp = "PH"),
             prior(student_t(3,0,1), class = "sd", resp = "Ell"),
             prior(student_t(3,0,1), class = "sd", resp = "PH"),
             prior(student_t(3,0,1), class = "sigma", resp = "Ell"),
             prior(student_t(3,0,1), class = "sigma", resp = "PH"),
             prior(student_t(3,0,5), class = "meanme", resp = "Ell"),
             prior(student_t(3,0,5), class = "meanme", resp = "PH"),
             prior(student_t(3,0,5), class = "sdme", resp = "Ell"),
             prior(student_t(3,0,5), class = "sdme", resp = "PH"),
             prior(lkj(5), class = "cor"))
lin_mod_hab <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH)*Management + Ndep*Management + 
                            (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                         bf(PH | mi(PH_SE)  ~ fieldseason_rain + Management*Sdep + 
                              Management*Ndep*Year1_pH + 
                              (1|p|SQUARE) +
                              ar(time = YRnm, gr = REP_ID)) + 
                         set_rescor(FALSE),
                       data = mod_data, prior = mod_pr,
                       cores = 4, iter = 5000, seed = 6856231456,
                       file = "Outputs/Models/Difference/Test_linmodhab_130421")
summary(lin_mod_hab)
plot(conditional_effects(lin_mod_hab))


# simplified linear
lin_mod_hab2 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH)*Management + 
                        (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                     bf(PH | mi(PH_SE)  ~ fieldseason_rain + Management*Sdep + 
                          Management*Ndep + 
                          (1|p|SQUARE) +
                          ar(time = YRnm, gr = REP_ID)) + 
                     set_rescor(FALSE),
                   data = mod_data, prior = mod_pr,
                   cores = 4, iter = 5000, seed = 6856231456,
                   file = "Outputs/Models/Difference/Test_linmodhab2_130421")
summary(lin_mod_hab2)
plot(conditional_effects(lin_mod_hab2))


lin_mod_hab3 <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH)*Management + 
                         (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) +
                      bf(PH | mi(PH_SE)  ~ fieldseason_rain + Management*Sdep + 
                           Management*Ndep + Year1_pH +
                           (1|p|SQUARE) +
                           ar(time = YRnm, gr = REP_ID)) + 
                      set_rescor(FALSE),
                    data = mod_data, prior = mod_pr,
                    cores = 4, iter = 5000, seed = 6856231456,
                    file = "Outputs/Models/Difference/Test_linmodhab3_130421")
summary(lin_mod_hab3)
plot(conditional_effects(lin_mod_hab3))


# compare predicted effects of Sdep and Ndep upon pH with/out splines or initial
# pH effects
splplots <- plot(conditional_effects(full_mod_hab_v4), points = TRUE, plot = FALSE,
                 point_args = list(alpha = 0.3))
linplots <- plot(conditional_effects(lin_mod_hab2), points = TRUE, plot = FALSE,
                 point_args = list(alpha = 0.3))
splplots[["PH.PH_Ndep:Management"]] + linplots[["PH.PH_Ndep:Management"]]
splplots[["PH.PH_Sdep:Management"]] + linplots[["PH.PH_Sdep:Management"]]

splplots2 <- plot(conditional_effects(full_mod_hab_v3), points = TRUE, plot = FALSE,
                  point_args = list(alpha = 0.3))
linplots2 <- plot(conditional_effects(lin_mod_hab3), points = TRUE, plot = FALSE,
                  point_args = list(alpha = 0.3))
splplots2[["PH.PH_Ndep:Management"]] + linplots2[["PH.PH_Ndep:Management"]]
splplots2[["PH.PH_Sdep:Management"]] + linplots2[["PH.PH_Sdep:Management"]]

p1 <- (splplots[["PH.PH_Ndep:Management"]] + linplots[["PH.PH_Ndep:Management"]] & 
  ggtitle("No initial pH effect")) 
p2 <- (splplots2[["PH.PH_Ndep:Management"]] + linplots2[["PH.PH_Ndep:Management"]] &
  ggtitle("Initial pH effect"))

p1/p2

p1 <- (splplots[["PH.PH_Sdep:Management"]] + linplots[["PH.PH_Sdep:Management"]] & 
         ggtitle("No initial pH effect")) 
p2 <- (splplots2[["PH.PH_Sdep:Management"]] + linplots2[["PH.PH_Sdep:Management"]] &
         ggtitle("Initial pH effect"))

p1/p2

# ** No measurement error ####
# ~~~ 200m2 weighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         N = NC_RATIO) %>%
  mutate(Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         N = as.numeric(scale(N)), 
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()

# Multivariate with N as a 2D spline
# Habitat interaction
get_prior(bf(Ell  ~ Management*PH + 
               Management*Ndep + Management*Sdep +
               (1|p|SQUARE) +
               ar(time = YRnm, gr = REP_ID)) +
            bf(PH  ~ Management*Sdep + fieldseason_rain + 
                 Management*Ndep + Year1_pH +
                 (1|p|SQUARE) +
                 ar(time = YRnm, gr = REP_ID)) + 
            set_rescor(FALSE), data = mod_data)


mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(normal(2,0.5), class = "ar", resp = "Ell"),
            prior(normal(2,0.5), class = "ar", resp = "PH"),
            prior(student_t(3,0,1), class = "sd", resp = "Ell"),
            prior(student_t(3,0,1), class = "sd", resp = "PH"),
            prior(student_t(3,0,1), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,1), class = "sigma", resp = "PH"),
            prior(lkj(5), class = "cor"))


# prior simulation
prior_mod <- brm(bf(Ell  ~ Management*PH + 
                      Management*Ndep + Management*Sdep +
                      (1|p|SQUARE) +
                      ar(time = YRnm, gr = REP_ID)) +
                   bf(PH  ~ Management*Sdep + fieldseason_rain + 
                        Management*Ndep + Year1_pH +
                        (1|p|SQUARE) +
                        ar(time = YRnm, gr = REP_ID)) + 
                   set_rescor(FALSE), data = mod_data, prior = mod_pr,
                 sample_prior = "only", save_pars = save_pars(all = TRUE, latent = TRUE), 
                 cores = 4, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50, resp = "Ell")
pp_check(prior_mod, nsamples = 50, resp = "PH")
plot(conditional_effects(prior_mod))

# run model - full weighted Ell R
mod_whw <- brm(bf(Ell  ~ Management*PH + 
                    (1|p|SQUARE) +
                    ar(time = YRnm, gr = REP_ID)) +
                 bf(PH  ~ Management*Sdep + fieldseason_rain + 
                      Management*Ndep + Year1_pH +
                      (1|p|SQUARE) +
                      ar(time = YRnm, gr = REP_ID)) + 
                 set_rescor(FALSE), data = mod_data, prior = mod_pr,
               save_pars = save_pars(all = TRUE, latent = TRUE), 
               iter = 20000, thin = 4, seed = 39099630,
               file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_WHW_PH_HAB",
               cores = 4)
summary(mod_whw)
plot(mod_whw)

mod_n_whw <- brm(bf(Ell  ~ Management*PH + 
                      Management*Ndep + 
                      (1|p|SQUARE) +
                      ar(time = YRnm, gr = REP_ID)) +
                   bf(PH  ~ Management*Sdep + fieldseason_rain + 
                        Management*Ndep + Year1_pH +
                        (1|p|SQUARE) +
                        ar(time = YRnm, gr = REP_ID)) + 
                   set_rescor(FALSE), data = mod_data, prior = mod_pr,
                 save_pars = save_pars(all = TRUE, latent = TRUE), 
                 iter = 20000, thin = 4, seed = 39099630,
                 file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_N_WHW_PH_HAB",
                 cores = 4)
summary(mod_n_whw)
plot(mod_n_whw)

mod_ns_whw <- brm(bf(Ell  ~ Management*PH + 
                       Management*Ndep + Management*Sdep +
                       (1|p|SQUARE) +
                       ar(time = YRnm, gr = REP_ID)) +
                    bf(PH  ~ Management*Sdep + fieldseason_rain + 
                         Management*Ndep + Year1_pH +
                         (1|p|SQUARE) +
                         ar(time = YRnm, gr = REP_ID)) + 
                    set_rescor(FALSE), data = mod_data, prior = mod_pr,
                  save_pars = save_pars(all = TRUE, latent = TRUE), 
                  iter = 20000, thin = 4, seed = 39099630,
                  file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_NS_WHW_PH_HAB",
                  cores = 4)
summary(mod_ns_whw)
plot(mod_ns_whw)

# ~~~ 200m2 unweighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         PH_SE = PH_DIW_SE_NORM,
         ELL_SE = ELL_WH_UW_SE_NORM) %>%
  mutate(Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()

# run model - full unweighted Ell R
mod_whuw <- update(mod_whw, newdata = mod_data,
                   cores = 4, iter = 20000, thin = 4, seed = 39099630,
                   file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_WHUW_PH_HAB")
summary(mod_whuw)
plot(mod_whuw)

mod_n_whuw <- update(mod_n_whw, newdata = mod_data,
                     cores = 4, iter = 20000, thin = 4, seed = 39099630,
                     file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_N_WHUW_PH_HAB")
summary(mod_n_whuw)
plot(mod_n_whuw)

mod_ns_whuw <- update(mod_ns_whw, newdata = mod_data,
                      cores = 4, iter = 20000, thin = 4, seed = 39099630,
                      file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_NS_WHUW_PH_HAB")
summary(mod_ns_whuw)
plot(mod_ns_whuw)

# ~~~ 4m2 weighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         PH_SE = PH_DIW_SE_NORM,
         ELL_SE = ELL_SM_W_SE_NORM) %>%
  mutate(Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()

# Multivariate with N as a 2D spline
# Habitat interaction
# run model - small weighted Ell R
mod_smw <- update(mod_whw, newdata = mod_data,
                  cores = 4, iter = 20000, thin = 4, seed = 39099630,
                  file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_SMW_PH_HAB")
summary(mod_smw)
plot(mod_smw)

mod_n_smw <- update(mod_n_whw, newdata = mod_data,
                    cores = 4,  iter = 20000, thin = 4, seed = 39099630,
                    file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_N_SMW_PH_HAB")
summary(mod_n_smw)
plot(mod_n_smw)

mod_ns_smw <- update(mod_ns_whw, newdata = mod_data,
                     cores = 4,  iter = 20000, thin = 4, seed = 39099630,
                     file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_NS_SMW_PH_HAB")
summary(mod_ns_smw)
plot(mod_ns_smw)


# ~~~ 4m2 unweighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         PH_SE = PH_DIW_SE_NORM,
         ELL_SE = ELL_SM_UW_SE_NORM) %>%
  mutate(Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()

# Multivariate with N as a 2D spline
# Habitat interaction
# run model - small unweighted Ell R
mod_smuw <- update(mod_whw, newdata = mod_data,
                   cores = 4, iter = 20000, thin = 4, seed = 39099630,
                   file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_SMUW_PH_HAB")
summary(mod_smuw)
plot(mod_smuw)

mod_n_smuw <- update(mod_n_whw, newdata = mod_data,
                     cores = 4, iter = 20000, thin = 4, seed = 39099630,
                     file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_N_SMUW_PH_HAB")
summary(mod_n_smuw)
plot(mod_n_smuw)

mod_ns_smuw <- update(mod_ns_whw, newdata = mod_data,
                      cores = 4, iter = 20000, thin = 4, seed = 39099630,
                      file = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/ELL_NS_SMUW_PH_HAB")
summary(mod_ns_smuw)
plot(mod_ns_smuw)

# test model with different time periods
mod_pr_time <- c(prior(normal(0,0.5), class = "b"),
                 prior(student_t(3, 0, 2.5), class = "sds"),
                 prior(normal(0,0.25), class = "Intercept"),
                 prior(gamma(4,1), class = "nu"),
                 prior(student_t(3,0,0.5), class = "sd"),
                 prior(student_t(3,0,0.5), class = "sigma"))

test_mod_yr1 <- brm(Ell  ~ Management*PH + s(Year1_pH, N, Management) +
                      (1|SQUARE),family = "student",
                    data = filter(mod_data, YRnm == 1),
                    prior = mod_pr_time,
                    save_pars = save_pars(all = TRUE),
                    cores = 4, iter = 5000,
                    control = list(adapt_delta = 0.99))
summary(test_mod_yr1)

test_mod_yr2 <- update(test_mod_yr1,
                       newdata = filter(mod_data, YRnm == 2),
                       save_pars = save_pars(all = TRUE),
                       cores = 4, iter = 5000,
                       control = list(adapt_delta = 0.99))
summary(test_mod_yr2)

test_mod_yr3 <- update(test_mod_yr1,
                       newdata = filter(mod_data, YRnm == 3),
                       save_pars = save_pars(all = TRUE),
                       cores = 4, iter = 5000,
                       control = list(adapt_delta = 0.99))
summary(test_mod_yr3)
plot(conditional_effects(test_mod_yr3,
                         conditions = data.frame(cond__ = c("Management = 1","Management = 0"))))

nd <-
  tibble(Year1_pH = seq(from = -2.25, to = 3.6, length.out = 30) %>%
           rep(., times = 6),
         N = rep(c(-1,0,1), each = 30) %>%
           rep(., times = 2),
         Management = rep(0:1, each = 90),
         PH = 0,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(test_mod_yr1, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Time_period = "7898"),
  fitted(test_mod_yr3, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Time_period = "9807"),
  fitted(test_mod_yr3, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Time_period = "0719")
)) %>%
  mutate(Time_period = forcats::fct_relevel(
    as.factor(Time_period),"7898","9807","0719"))

plot_dat <- ELL_pH %>%
  select(ELL = SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Year1_pH, N = NC_RATIO) %>%
  filter(!is.na(Management)) %>%
  mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"))

plot_dat %>%
  mutate(N = cut(scale(N), c(-10,-0.5,0.5,10), labels = c("-1","0","1")),
         Year1_pH = scale(Year1_pH),
         Time_period = as.factor(Time_period)) %>%
  mutate(Time_period = forcats::fct_relevel(Time_period,"7898","9807","0719")) %>%
  filter(!is.na(N)) %>%
  ggplot() +
  geom_point(aes(x = Year1_pH, y = ELL, colour = N),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity",
              alpha = 1/3, size = 1/2) +
  facet_grid(Time_period~Management) +
  labs(x = "Initial pH", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Initial pH and N effect on Ellenberg R over different time periods with data.png",
       path = "Outputs/Models/Difference/Univariate_nomeaserror",
       width = 16, height = 18, units = "cm")

nd <-
  tibble(PH = seq(from = -3, to = 3.5, length.out = 30) %>%
           rep(., times = 2),
         N = 0,
         Management = rep(0:1, each = 30),
         Year1_pH = 0,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(test_mod_yr1, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Time_period = "7898"),
  fitted(test_mod_yr3, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Time_period = "9807"),
  fitted(test_mod_yr3, newdata = nd, re_formula = NA) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Time_period = "0719")
)) %>%
  mutate(Time_period = forcats::fct_relevel(
    as.factor(Time_period),"7898","9807","0719"))

plot_dat %>%
  mutate(Time_period = as.factor(Time_period)) %>%
  mutate(Time_period = forcats::fct_relevel(Time_period,"7898","9807","0719")) %>%
  filter(!is.na(N)) %>%
  ggplot() +
  geom_point(aes(x = PH, y = ELL, colour = Management),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity",
              alpha = 1/3, size = 1/2) +
  facet_grid(~Time_period) +
  labs(x = "Initial pH", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("pH change on Ellenberg R over different time periods with data.png",
       path = "Outputs/Models/Difference/Univariate_nomeaserror",
       width = 20, height = 12, units = "cm")


# kfold ####
options(future.globals.maxSize = 10e8)
library(future)
plan(multiprocess)

mods_list <- list.files(path = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/",
                        pattern = ".rds$")

for (i in 1:length(mods_list)){
  mod <- readRDS(mods_list[i])
  
  filename <- gsub(".rds","",mods_list[i])
  
  mod <- update(mod, cores = 4, iter = 20000, thin = 4,
                seed = 39099630)
  
  mod <- add_criterion(mod, "kfold", group = "SQUARE",
                       folds = "grouped",
                       file = filename, overwrite = TRUE)
  
}

loo_compare(mod_whw, mod_whw_ns, criterion = "kfold")
# elpd_diff se_diff
# mod_whw       0.0       0.0 
# mod_whw_ns -185.5      32.1 

loo_compare(mod_whuw, mod_whuw_ns, criterion = "kfold")
# elpd_diff se_diff
# mod_whuw       0.0       0.0 
# mod_whuw_ns -250.6      52.4 

loo_compare(mod_smw, mod_smw_ns, criterion = "kfold")
# elpd_diff se_diff
# mod_smw       0.0       0.0 
# mod_smw_ns -344.7      50.3 

loo_compare(mod_smuw, mod_smuw_ns, criterion = "kfold")
# elpd_diff se_diff
# mod_smuw     0.0       0.0   
# mod_smuw_ns -7.1      11.5   


# Summary plots ####
# Plot for comparing initial pH effects on pH change
nd <- 
  tibble(Year1_pH = seq(from = -2.25, to = 3.6, length.out = 30) %>% 
           rep(., times = 2),
         Management = rep(c("Low","High"), each = 30),
         Sdep = 0,
         Ndep = 0,
         fieldseason_rain = 0,
         PH_SE = 0.3,
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted small")
)) %>%
  mutate(Year1_pH = 1.2978*Year1_pH + 5.598)

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Year1_pH, Ndep) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))

plot_dat %>%
  ggplot() +
  geom_point(aes(x = Year1_pH, y = PH),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = filter(f, Management == "Low intensity"),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity",
              alpha = 1/3, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Initial pH", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Initial pH effect on pH with data measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models", 
       width = 15, height = 12, units = "cm")

# Plot for comparing Sdep effects on pH change
nd <- 
  tibble(Sdep = seq(from = -5.5, to = 1.75, length.out = 30) %>% 
           rep(., times = 2),
         fieldseason_rain = 0,
         Management = rep(c("Low","High"), each = 30),
         Year1_pH = 0,
         Ndep = 0,
         PH_SE = 0.3,
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted small")
)) %>%
  mutate(
    Sdep = 3.818393*Sdep - 4.769306
  )

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Year1_pH, N = NC_RATIO,
         Sdep, fieldseason_rain) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))


plot_dat %>%
  ggplot() +
  geom_point(aes(x = Sdep, y = PH, colour = Management),
             alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity", 
              alpha = 1/2, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("Sdep effect on pH with data measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/", 
       width = 16, height = 10, units = "cm")

# mean and prediction error
f2 <- do.call(rbind, list(
  posterior_predict(mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    summarise(across(.fns = list(Estimate = mean,
                                 Q2.5 = ~quantile(.x, probs = 0.025),
                                 Q97.5 = ~quantile(.x, probs = 0.975)))) %>%
    pivot_longer(everything(), names_to = c("ID","Metric"),
                 names_sep = "_") %>%
    pivot_wider(names_from = "Metric", values_from = "value") %>%
    bind_cols(nd) %>%
    mutate(Response = "Weighted full"),
  posterior_predict(mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    summarise(across(.fns = list(Estimate = mean,
                                 Q2.5 = ~quantile(.x, probs = 0.025),
                                 Q97.5 = ~quantile(.x, probs = 0.975)))) %>%
    pivot_longer(everything(), names_to = c("ID","Metric"),
                 names_sep = "_") %>%
    pivot_wider(names_from = "Metric", values_from = "value") %>%
    bind_cols(nd) %>%
    mutate(Response = "Unweighted full"),
  posterior_predict(mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    summarise(across(.fns = list(Estimate = mean,
                                 Q2.5 = ~quantile(.x, probs = 0.025),
                                 Q97.5 = ~quantile(.x, probs = 0.975)))) %>%
    pivot_longer(everything(), names_to = c("ID","Metric"),
                 names_sep = "_") %>%
    pivot_wider(names_from = "Metric", values_from = "value") %>%
    bind_cols(nd) %>%
    mutate(Response = "Weighted small"),
  posterior_predict(mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    summarise(across(.fns = list(Estimate = mean,
                                 Q2.5 = ~quantile(.x, probs = 0.025),
                                 Q97.5 = ~quantile(.x, probs = 0.975)))) %>%
    pivot_longer(everything(), names_to = c("ID","Metric"),
                 names_sep = "_") %>%
    pivot_wider(names_from = "Metric", values_from = "value") %>%
    bind_cols(nd) %>%
    mutate(Response = "Unweighted small")
)) %>%
  mutate(
    Management = ifelse(Management == "High", "High intensity", "Low intensity"),
    Sdep = 3.818393*Sdep - 4.769306
  )
plot_dat %>%
  ggplot() +
  geom_point(aes(x = Sdep, y = PH, colour = Management),
             alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_ribbon(data = f2,
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, 
                  group = Management),
              stat = "identity", 
              alpha = 1/4) +
  geom_smooth(data = f,
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity", 
              alpha = 1/2, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0))

# Plot for comparing Ndep effects on pH change
nd <- 
  tibble(Ndep = seq(from = -1.8, to = 3.15, length.out = 30) %>% 
           rep(., times = 2),
         fieldseason_rain = 0,
         Management = rep(c("Low","High"), each = 30),
         Year1_pH = 0,
         Sdep = 0,
         PH_SE = 0.3,
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted small")
)) %>%
  mutate(
    Ndep = 70.08984*Ndep + 138.8274
  )

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Year1_pH, Ndep,
         Sdep, fieldseason_rain) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))


plot_dat %>%
  ggplot() +
  geom_point(aes(x = Ndep, y = PH, colour = Management),
             alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity", 
              alpha = 1/2, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "N deposition", y = "pH change") +
  scale_x_continuous(expand = expansion(),limits = c(0,390)) 
ggsave("Ndep effect on pH with data measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/", 
       width = 16, height = 10, units = "cm")

# Plot for comparing pH change on Ellenberg R change
nd <- 
  tibble(PH = seq(from = -3, to = 3.5, length.out = 40) %>% 
           rep(., times = 2),
         Management = rep(c("High","Low"), each = 40),
         Year1_pH = 0,
         N = 0,
         ELL_SE = 0.4,
         YRnm = 1,
         REP_ID = 1:80)

f <- do.call(rbind, list(
  fitted(mod_whw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(mod_whuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(mod_smw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(mod_smuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Year1_pH, N = NC_RATIO) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_grid(~Response) +
  labs(x = "pH change", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("pH change effect on Ellenberg R measerror.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/", 
       width = 20, height = 10, units = "cm")

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = PH, y = value, colour = Management),
             alpha = 0.5) +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(~Response) +
  labs(x = "pH change", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("pH change effect on Ellenberg R with data measerror.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models", 
       width = 20, height = 10, units = "cm")

# Parameter estimates
mcmc_plot(mod_whw, type = "areas_ridges") + 
  ggtitle(expression("Cover weighted Ellenberg R 200m"^2)) + 
  mcmc_plot(mod_whuw, type = "areas_ridges") + theme(axis.text.y = element_blank())+
  ggtitle(expression("Unweighted Ellenberg R 200m"^2)) +
  mcmc_plot(mod_smw, type = "areas_ridges") + 
  ggtitle(expression("Cover weighted Ellenberg R 4m"^2)) +
  mcmc_plot(mod_smuw, type = "areas_ridges") + theme(axis.text.y = element_blank()) +
  ggtitle(expression("Unweighted Ellenberg R 4m"^2)) 

ggsave("Parameter estimates bayesplot.png", 
       path = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/",
       width = 24, height = 14, units = "cm", scale = 1.2)


plot_pars <- c(parnames(mod_whw)[1:21],"ar_PH","ar_Ell")
param_summaries <- do.call(rbind, list(
  posterior_summary(mod_whw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "WH_W"),
  posterior_summary(mod_whuw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "WH_UW"),
  posterior_summary(mod_smw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "SM_W"),
  posterior_summary(mod_smuw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "SM_UW")
)) 
write.csv(param_summaries, 
          "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/Parameter_summary_table.csv",
          row.names = FALSE)

param_summaries <- param_summaries %>% 
  mutate(across(Estimate:Q97.5, round, 3)) %>%
  mutate(CI = paste0(Q2.5, " - ", Q97.5),
         Estimate = as.character(Estimate)) %>%
  select(-Est.Error, -Q2.5, -Q97.5) %>%
  pivot_longer(c(Estimate,CI)) %>%
  pivot_wider(names_from = c(Model, name), names_sep = "__") %>%
  mutate(across(ends_with("Estimate"),as.numeric))
writexl::write_xlsx(param_summaries, 
                    "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/Parameter_summary_nicetable.xlsx")


# NS as Ellenberg predictors
mcmc_plot(mod_ns_whw, type = "areas_ridges") + 
  ggtitle(expression("Cover weighted Ellenberg R 200m"^2)) + 
  mcmc_plot(mod_ns_whuw, type = "areas_ridges") + theme(axis.text.y = element_blank())+
  ggtitle(expression("Unweighted Ellenberg R 200m"^2)) +
  mcmc_plot(mod_ns_smw, type = "areas_ridges") + 
  ggtitle(expression("Cover weighted Ellenberg R 4m"^2)) +
  mcmc_plot(mod_ns_smuw, type = "areas_ridges") + theme(axis.text.y = element_blank()) +
  ggtitle(expression("Unweighted Ellenberg R 4m"^2)) 

ggsave("Parameter estimates Ell_NS models bayesplot.png", 
       path = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/",
       width = 24, height = 14, units = "cm", scale = 1.3)


plot_pars <- c(parnames(mod_ns_whw)[1:25],"ar_PH","ar_Ell")
param_summaries <- do.call(rbind, list(
  posterior_summary(mod_ns_whw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "WH_W"),
  posterior_summary(mod_ns_whuw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "WH_UW"),
  posterior_summary(mod_ns_smw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "SM_W"),
  posterior_summary(mod_ns_smuw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "SM_UW")
)) %>% 
  mutate(across(Estimate:Q97.5, round, 3)) %>%
  mutate(CI = paste0(Q2.5, " - ", Q97.5),
         Estimate = as.character(Estimate)) %>%
  select(-Est.Error, -Q2.5, -Q97.5) %>%
  pivot_longer(c(Estimate,CI)) %>%
  pivot_wider(names_from = c(Model, name), names_sep = "__") %>%
  mutate(across(ends_with("Estimate"),as.numeric))
writexl::write_xlsx(param_summaries, 
                    "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models/Parameter_summary_nsmodel_nicetable.xlsx")

plot_pars <- c(parnames(mod_whw)[1:21],"ar_PH","ar_Ell")
post <- do.call(rbind, list(
  posterior_samples(mod_whw, pars = plot_pars) %>% 
    mutate(Model = "WH_W"),
  posterior_samples(mod_whuw, pars = plot_pars) %>% 
    mutate(Model = "WH_UW"),
  posterior_samples(mod_smw, pars = plot_pars) %>% 
    mutate(Model = "SM_W"),
  posterior_samples(mod_smuw, pars = plot_pars) %>% 
    mutate(Model = "SM_UW")
)) 
post %>%
  transmute(Ell_pH_High    = b_Ell_PH + `b_Ell_Management:PH`,
            Ell_pH_Low = b_Ell_PH,
            pH_Sdep_High = b_PH_Sdep + `b_PH_Management:Sdep`,
            pH_Sdep_Low = b_PH_Sdep,
            pH_Ndep_High = b_PH_Ndep + `b_PH_Management:Ndep`,
            pH_Ndep_Low = b_PH_Ndep,
            # pH_InitialpH_NA = b_PH_Year1_pH,
            # pH_Rainfall_NA = b_PH_fieldseason_rain,
            Model = Model) %>%
  tidyr::pivot_longer(contains("_"), "key", "value") %>%
  tidyr::separate(key, c("Response","Predictor","Management"),
                  sep = "_", remove = FALSE) %>%
  mutate(group = paste0(Model,key),
         Response = paste("Response:",Response),
         Predictor = paste("Predictor:",Predictor),
         Model = recode(Model,
                        "SM_UW" = "Small unweighted",
                        "SM_W" = "Small weighted",
                        "WH_W" = "Whole weighted",
                        "WH_UW" = "Whole unweighted")) %>%
  ggplot(aes(x = value, group = group, 
             color = Management, fill = Management)) +
  geom_density(alpha = 1/4) +
  facet_grid(Model~Response + Predictor) +
  scale_x_continuous("Parameter estimate", expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0,0)) +
  scale_colour_manual(values = mang_cols, aesthetics = c("colour","fill")) 
ggsave("Parameter estimates of pH and Ndep and Sdep by management.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror/Linear_models",
       width = 15, height = 15, units = "cm")


# ** Measurement error ####
# Multivariate with Ndep effects upon pH only
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         ELL_SE = ELL_WH_W_SE_NORM, 
         PH_SE = PH_DIW_SE_NORM,
         N = NC_RATIO) %>%
  mutate(Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)),
         N = as.numeric(scale(N)),
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()


# Habitat interaction
get_prior(bf(Ell | mi(ELL_SE) ~ Management*mi(PH) +  
               (1|p|SQUARE) +
               ar(time = YRnm, gr = REP_ID)) +
            bf(PH | mi(PH_SE)  ~ Management*Sdep + fieldseason_rain + 
                 Management*Ndep + Year1_pH + 
                 (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) + 
            set_rescor(FALSE), data = mod_data)


mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(normal(2,0.5), class = "ar", resp = "Ell"),
            prior(normal(2,0.5), class = "ar", resp = "PH"),
            prior(student_t(3,0,1), class = "sd", resp = "Ell"),
            prior(student_t(3,0,1), class = "sd", resp = "PH"),
            prior(student_t(3,0,1), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,1), class = "sigma", resp = "PH"),
            prior(student_t(3,0,5), class = "meanme", resp = "Ell"),
            prior(student_t(3,0,5), class = "meanme", resp = "PH"),
            prior(student_t(3,0,5), class = "sdme", resp = "Ell"),
            prior(student_t(3,0,5), class = "sdme", resp = "PH"),
            prior(lkj(5), class = "cor"))

# prior simulation
prior_mod <- brm(bf(Ell | mi(ELL_SE) ~ Management*mi(PH) +  
                      (1|p|SQUARE) +
                      ar(time = YRnm, gr = REP_ID)) +
                   bf(PH | mi(PH_SE)  ~ Management*Sdep + fieldseason_rain + 
                        Management*Ndep + Year1_pH + 
                        (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) + 
                   set_rescor(FALSE), data = mod_data, prior = mod_pr,
                 sample_prior = "only", 
                 save_pars = save_pars(all = TRUE, latent = TRUE), 
                 cores = 4, iter = 5000)
summary(prior_mod)
plot(prior_mod)
pp_check(prior_mod, nsamples = 50, resp = "Ell")
pp_check(prior_mod, nsamples = 50, resp = "PH")
plot(conditional_effects(prior_mod))

# simulate data according to model and check if works 
set.seed(151020)
sq_devs <- rnorm(100,0,0.2)
sim_data <- data.frame(SQUARE = gl(100,5)) %>%
  mutate(REP_ID = paste0(SQUARE,"X",1:500),    
         Improved = rbinom(500,1,0.007*as.numeric(SQUARE))) %>%
  mutate(pH = rnorm(500,mean = 4.5 + Improved + rep(rnorm(100,sq_devs,0.2),each = 5), 1)) %>%
  mutate(Ell = rnorm(500, mean = pH + rep(rnorm(100,sq_devs,0.3),each = 5), 1)) %>%
  mutate(Sdep = rep(rnorm(100,0,1),each = 5)) %>%
  mutate(Ndep = rep(rnorm(100,-Sdep,1),each = 5)) %>%
  mutate(rain_year2 = rep(rnorm(100,0,1),each = 5),
         pHSE_year1 = 0.2,
         EllSE_year1 = 0.15 + Improved*0.05) %>%
  # mutate(test = 0.5*Improved+(1-Improved)*0.5*Sdep + 0.3*rain)
  mutate(pH_diffYear2 = rstudent_t(500, 4,
                                   0.5*Improved+(1-Improved)*0.5*Sdep + 
                                     0.3*rain_year2,pHSE_year1))%>%
  mutate(Ell_diffYear2 = rstudent_t(500, 4,
                                    pH_diffYear2*(1-Improved) + 
                                      Improved*0.5*pH_diffYear2,
                                    EllSE_year1)) %>%
  mutate(pH_year2 = pH + pH_diffYear2,
         pHSE_year2 = 0.15,
         Ell_year2 = Ell + Ell_diffYear2,
         EllSE_year2 = 0.2 + Improved*0.05,
         rain_year3 = rep(rnorm(100,0,1),each = 5)) %>%
  mutate(pH_diffYear3 = rstudent_t(500, 4,
                                   0.5*Improved+(1-Improved)*0.5*Sdep + 
                                     0.3*rain_year3, pHSE_year2)) %>%
  mutate(Ell_diffYear3 = rstudent_t(500, 4,
                                    pH_diffYear3*(1-Improved) + 
                                      Improved*0.5*pH_diffYear3,
                                    EllSE_year2)) %>%
  select(SQUARE, REP_ID, Sdep, Ndep, Improved,
         rain_year12 = rain_year2,rain_year23 = rain_year3,
         PH_SE_year12 = pHSE_year1, PH_SE_year23 = pHSE_year2,
         ELL_SE_year12 = EllSE_year1, ELL_SE_year23 = EllSE_year2,
         PH_year12 = pH_diffYear2, PH_year23 = pH_diffYear3,
         Year1_pH_year12 = pH, Year1_pH_year23 = pH_year2,
         Ell_year12 = Ell_diffYear2,Ell_year23 = Ell_diffYear3) %>%
  pivot_longer(contains("_year"), names_to = c("Variable","Time_period"),
               names_sep = "_year") %>%
  pivot_wider(names_from = Variable, values_from = value) %>%
  mutate(YRnm = as.numeric(as.factor(Time_period)))
psych::multi.hist(select_if(sim_data, is.numeric))
psych::pairs.panels(select_if(sim_data, is.numeric), rug = FALSE)


sim_mod <- brm(bf(Ell | mi(ELL_SE) ~ Improved*mi(PH) +  
                    (1|p|SQUARE) +
                    ar(time = YRnm, gr = REP_ID)) +
                 bf(PH | mi(PH_SE)  ~ Improved*Sdep + rain + (1|p|SQUARE) +
                      ar(time = YRnm, gr = REP_ID)) + 
                 set_rescor(FALSE), data = sim_data, prior = mod_pr,
               save_pars = save_pars(all = TRUE, latent = TRUE), 
               cores = 4, iter = 4000)
summary(sim_mod)
pp_check(sim_mod, nsamples = 50, resp = "Ell")
pp_check(sim_mod, nsamples = 50, resp = "PH")
pp_check(sim_mod, nsamples = 50, resp = "N")

plot(conditional_effects(sim_mod))


# ~~~ 200m2 weighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         ELL_SE = ELL_WH_W_SE_NORM, 
         PH_SE = PH_DIW_SE_NORM) %>%
  mutate(Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()


# run model - full weighted Ell R
mod_whw <- brm(bf(Ell | mi(ELL_SE)  ~ Management*mi(PH) +
                         (1|p|SQUARE) +
                         ar(time = YRnm, gr = REP_ID)) +
                      bf(PH | mi(PH_SE)  ~ Management*Sdep + fieldseason_rain  + 
                           Management*Ndep + Year1_pH +
                           (1|p|SQUARE) +
                           ar(time = YRnm, gr = REP_ID)) + 
                      set_rescor(FALSE), data = mod_data, prior = mod_pr,
                    save_pars = save_pars(all = TRUE, latent = TRUE), 
                    file = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/ELL_WHW_PH_HAB",
                    cores = 4, iter = 20000, thin = 4)
summary(mod_whw)
plot(mod_whw)
pp_check(mod_whw, nsamples = 50, resp = "Ell")
pp_check(mod_whw, nsamples = 50, resp = "PH")

# check residuals against habitat
mod_whw_resid <- resid(mod_whw)
dimnames(mod_whw_resid)
BH_plot <- left_join(mod_data,
                     ELL_pH %>%
                       filter(!is.na(Management)) %>%
                       mutate(YRnm = as.integer(as.factor(Year))) %>%
                       select(REP_ID, YRnm, BH, BH_DESC) ) %>%
  mutate(BH_DESC = recode(BH_DESC,
                          "Broadleaved Mixed and Yew Woodland" = "Broadleaved Woodland"))
all.equal(BH_plot$REP_ID, mod_data$REP_ID)
all.equal(BH_plot$YRnm, mod_data$YRnm)

BH_plot <- mutate(BH_plot,
                  Ell_resid = mod_whw_resid[,1,"Ell"],
                  PH_resid = mod_whw_resid[,1,"PH"]) %>%
  pivot_longer(ends_with("resid")) %>%
  mutate(name = sapply(strsplit(name, "_"),"[",1))

ggplot(BH_plot, aes(x = BH_DESC, y = value)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.5,
              colour = "dodgerblue3") +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Residuals") +
  facet_wrap(~name, ncol = 1) +
  theme(axis.text.x = element_text(angle = 60,hjust=1))
ggsave("Residuals by habitat whole weighted model.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models",
       width = 12, height = 18, units = "cm")

mod_whw_ns <- mod_whw <- brm(bf(Ell | mi(ELL_SE)  ~ Management*mi(PH) +
                                  Management*Sdep + Management*Ndep +
                                  (1|p|SQUARE) +
                                  ar(time = YRnm, gr = REP_ID)) +
                               bf(PH | mi(PH_SE)  ~ Management*Sdep + fieldseason_rain  + 
                                    Management*Ndep + Year1_pH +
                                    (1|p|SQUARE) +
                                    ar(time = YRnm, gr = REP_ID)) + 
                               set_rescor(FALSE), data = mod_data, prior = mod_pr,
                             save_pars = save_pars(all = TRUE, latent = TRUE), 
                             file = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/ELL_NS_WHW_PH_HAB",
                             cores = 4, iter = 20000, thin = 4)
summary(mod_whw_ns)
plot(mod_whw_ns)

# ~~~ 200m2 unweighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         PH_SE = PH_DIW_SE_NORM,
         ELL_SE = ELL_WH_UW_SE_NORM) %>%
  mutate(Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()

# Multivariate with N as a 2D spline
# Habitat interaction
# run model - full unweighted Ell R
mod_whuw <- update(mod_whw, newdata = mod_data,
                   cores = 4, iter = 20000, thin = 4,
                   file = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/ELL_WHUW_PH_HAB")

summary(mod_whuw)
plot(mod_whuw, ask = FALSE)
pp_check(mod_whuw, nsamples = 50, resp = "Ell")
pp_check(mod_whuw, nsamples = 50, resp = "PH")

# plot(conditional_effects(mod_whuw))
# check residuals against habitat
mod_resid <- resid(mod_whuw)
BH_plot <- left_join(mod_data,
                     ELL_pH %>%
                       filter(!is.na(Management)) %>%
                       mutate(YRnm = as.integer(as.factor(Year))) %>%
                       select(REP_ID, YRnm, BH, BH_DESC) ) %>%
  mutate(BH_DESC = recode(BH_DESC,
                          "Broadleaved Mixed and Yew Woodland" = "Broadleaved Woodland"))
BH_plot <- mutate(BH_plot,
                  Ell_resid = mod_resid[,1,"Ell"],
                  PH_resid = mod_resid[,1,"PH"]) %>%
  pivot_longer(ends_with("resid")) %>%
  mutate(name = sapply(strsplit(name, "_"),"[",1))

ggplot(BH_plot, aes(x = BH_DESC, y = value)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.5,
              colour = "dodgerblue3") +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Residuals") +
  facet_wrap(~name, ncol = 1) +
  theme(axis.text.x = element_text(angle = 60,hjust=1))
ggsave("Residuals by habitat whole unweighted model.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/",
       width = 12, height = 18, units = "cm")

mod_whuw_ns <- update(mod_whw_ns, newdata = mod_data,
                      cores = 4, iter = 20000, thin = 4,
                      file = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/ELL_NS_WHUW_PH_HAB")

summary(mod_whuw_ns)
plot(mod_whuw_ns, ask = FALSE)


# ~~~ 4m2 weighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         PH_SE = PH_DIW_SE_NORM,
         ELL_SE = ELL_SM_W_SE_NORM) %>%
  mutate(Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()

# Multivariate with N as a 2D spline
# Habitat interaction
# run model - small weighted Ell R
mod_smw <- update(mod_whw, newdata = mod_data,
                  cores = 4, iter = 20000, thin = 4,
                  file = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/ELL_SMW_PH_HAB")

summary(mod_smw)
plot(mod_smw)
pp_check(mod_smw, nsamples = 50, resp = "Ell")
pp_check(mod_smw, nsamples = 50, resp = "PH")

# plot(conditional_effects(mod_smw))
mod_resid <- resid(mod_smw)
BH_plot <- left_join(mod_data,
                     ELL_pH %>%
                       filter(!is.na(Management)) %>%
                       mutate(YRnm = as.integer(as.factor(Year))) %>%
                       select(REP_ID, YRnm, BH, BH_DESC) ) %>%
  mutate(BH_DESC = recode(BH_DESC,
                          "Broadleaved Mixed and Yew Woodland" = "Broadleaved Woodland"))
BH_plot <- mutate(BH_plot,
                  Ell_resid = mod_resid[,1,"Ell"],
                  PH_resid = mod_resid[,1,"PH"]) %>%
  pivot_longer(ends_with("resid")) %>%
  mutate(name = sapply(strsplit(name, "_"),"[",1))

ggplot(BH_plot, aes(x = BH_DESC, y = value)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.5,
              colour = "dodgerblue3") +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Residuals") +
  facet_wrap(~name, ncol = 1) +
  theme(axis.text.x = element_text(angle = 60,hjust=1))
ggsave("Residuals by habitat small weighted full model.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")

mod_smw_ns <- update(mod_whw_ns, newdata = mod_data,
                      cores = 4, iter = 20000, thin = 4,
                      file = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/ELL_NS_SMW_PH_HAB")

summary(mod_smw_ns)
plot(mod_smw_ns, ask = FALSE)

# ~~~ 4m2 unweighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year1_pH, PH, 
         PH_SE = PH_DIW_SE_NORM,
         ELL_SE = ELL_SM_UW_SE_NORM) %>%
  mutate(Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()

# Habitat interaction
# run model - small unweighted Ell R
mod_smuw <- update(mod_whw, newdata = mod_data,
                   cores = 4, iter = 20000, thin = 4,
                   file = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/ELL_SMUW_PH_HAB")

summary(mod_smuw)
plot(mod_smuw)
pp_check(mod_smuw, nsamples = 50, resp = "Ell")
pp_check(mod_smuw, nsamples = 50, resp = "PH")

# plot(conditional_effects(mod_smuw))

mod_resid <- resid(mod_smuw)
BH_plot <- left_join(mod_data,
                     ELL_pH %>%
                       filter(!is.na(Management)) %>%
                       mutate(YRnm = as.integer(as.factor(Year))) %>%
                       select(REP_ID, YRnm, BH, BH_DESC) ) %>%
  mutate(BH_DESC = recode(BH_DESC,
                          "Broadleaved Mixed and Yew Woodland" = "Broadleaved Woodland"))
BH_plot <- mutate(BH_plot,
                  Ell_resid = mod_resid[,1,"Ell"],
                  PH_resid = mod_resid[,1,"PH"]) %>%
  pivot_longer(ends_with("resid")) %>%
  mutate(name = sapply(strsplit(name, "_"),"[",1))

ggplot(BH_plot, aes(x = BH_DESC, y = value)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.5,
              colour = "dodgerblue3") +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Residuals") +
  facet_wrap(~name, ncol = 1) +
  theme(axis.text.x = element_text(angle = 60,hjust=1))
ggsave("Residuals by habitat small unweighted model.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/",
       width = 12, height = 18, units = "cm")


mod_smuw_ns <- update(mod_whw_ns, newdata = mod_data,
                      cores = 4, iter = 20000, thin = 4,
                      file = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/ELL_NS_SMUW_PH_HAB")
summary(mod_smuw_ns)
plot(mod_smuw_ns, ask = FALSE)

post <- posterior_samples(mod_smuw) 

post %>%
  transmute(Ell_pH_ManagementLow    = bsp_Ell_miPH + `bsp_Ell_miPH:ManagementLow`,
            Ell_pH_ManagementHigh = bsp_Ell_miPH,
            pH_Sdep_Low = b_PH_Sdep + `b_PH_ManagementLow:Sdep`,
            pH_Sdep_High = b_PH_Sdep,
            pH_Ndep_Low = b_PH_Ndep + `b_PH_ManagementLow:Ndep`,
            pH_Ndep_High = b_PH_Ndep) %>%
  tidyr::gather(key, value) %>%
  group_by(key) %>%
  summarise(mean = mean(value), 
            Q2.5 = quantile(value, probs = 0.025),
            Q97.5 = quantile(value, probs = 0.975))
# # A tibble: 6 x 4
# key                      mean    Q2.5  Q97.5
# <chr>                   <dbl>   <dbl>  <dbl>
# 1 Ell_pH_ManagementHigh 0.186    0.113  0.258 
# 2 Ell_pH_ManagementLow  0.146    0.0840 0.211 
# 3 pH_Ndep_High          0.214    0.0861 0.343 
# 4 pH_Ndep_Low           0.156    0.0735 0.238 
# 5 pH_Sdep_High          0.0394  -0.0386 0.114 
# 6 pH_Sdep_Low           0.00287 -0.0658 0.0708

post %>%
  transmute(Ell_pH_Low    = bsp_Ell_miPH + `bsp_Ell_miPH:ManagementLow`,
            Ell_pH_High = bsp_Ell_miPH,
            pH_Sdep_Low = b_PH_Sdep + `b_PH_ManagementLow:Sdep`,
            pH_Sdep_High = b_PH_Sdep,
            pH_Ndep_Low = b_PH_Ndep + `b_PH_ManagementLow:Ndep`,
            pH_Ndep_High = b_PH_Ndep,
            pH_InitialpH_NA = b_PH_Year1_pH,
            pH_Rainfall_NA = b_PH_fieldseason_rain) %>%
  tidyr::gather(key, value) %>%
  tidyr::separate(key, c("Response","Predictor","Management"),
                  sep = "_", remove = FALSE) %>%
  ggplot(aes(x = value, group = key, 
             color = Management, fill = Management)) +
  geom_density(alpha = 1/4) +
  facet_wrap(~Response + Predictor) +
  scale_x_continuous(expression(gamma), expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL)
  

# Kfold ####
options(future.globals.maxSize = 10e8)
library(future)
plan(multiprocess)

mods_list <- list.files(path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/",
                        pattern = ".rds$")

for (i in 1:length(mods_list)){
  mod <- readRDS(mods_list[i])
  
  filename <- gsub(".rds","",mods_list[i])
  
  mod <- update(mod, cores = 4, iter = 20000, thin = 4,
                seed = 39099630)
  
  mod <- add_criterion(mod, "kfold", group = "SQUARE",
                       folds = "grouped",
                       file = filename, overwrite = TRUE)
  
}


loo_compare(mod_whw, mod_whw_ns, criterion = "kfold")
# elpd_diff se_diff
# mod_whw      0.0       0.0  
# mod_whw_ns -21.5       8.7  


loo_compare(mod_whuw, mod_whuw_ns, criterion = "kfold")
# elpd_diff se_diff
# mod_whuw_ns   0.0       0.0  
# mod_whuw    -34.4      15.9  

loo_compare(mod_smw, mod_smw_ns, criterion = "kfold")
# elpd_diff se_diff
# mod_smw_ns   0.0       0.0  
# mod_smw    -20.4      10.7  

loo_compare(mod_smuw, mod_smuw_ns, criterion = "kfold")
# elpd_diff se_diff
# mod_smuw     0.0       0.0   
# mod_smuw_ns -7.1      11.5   

# Summary plots ####
# Plot for comparing initial pH effects on pH change
nd <- 
  tibble(Year1_pH = seq(from = -2.25, to = 3.6, length.out = 30) %>% 
           rep(., times = 2),
         Management = rep(c("Low","High"), each = 30),
         Sdep = 0,
         Ndep = 0,
         fieldseason_rain = 0,
         PH_SE = 0.3,
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted small")
)) %>%
  mutate(Year1_pH = 1.2978*Year1_pH + 5.598)

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Year1_pH, Ndep) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))

plot_dat %>%
  ggplot() +
  geom_point(aes(x = Year1_pH, y = PH),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = filter(f, Management == "Low intensity"),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity",
              alpha = 1/3, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Initial pH", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Initial pH effect on pH with data measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models", 
       width = 15, height = 12, units = "cm")

# Plot for comparing Sdep effects on pH change
nd <- 
  tibble(Sdep = seq(from = -5.5, to = 1.75, length.out = 30) %>% 
           rep(., times = 2),
         fieldseason_rain = 0,
         Management = rep(c("Low","High"), each = 30),
         Year1_pH = 0,
         Ndep = 0,
         PH_SE = 0.3,
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted small")
)) %>%
  mutate(
    Sdep = 3.818393*Sdep - 4.769306
  )

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Year1_pH, N = NC_RATIO,
         Sdep, fieldseason_rain) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))


plot_dat %>%
  ggplot() +
  geom_point(aes(x = Sdep, y = PH, colour = Management),
             alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity", 
              alpha = 1/2, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("Sdep effect on pH with data measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/", 
       width = 16, height = 10, units = "cm")

# mean and prediction error - small unweighted only
f2 <- posterior_predict(mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    summarise(across(.fns = list(Estimate = mean,
                                 Q2.5 = ~quantile(.x, probs = 0.025),
                                 Q97.5 = ~quantile(.x, probs = 0.975)))) %>%
    pivot_longer(everything(), names_to = c("ID","Metric"),
                 names_sep = "_") %>%
    pivot_wider(names_from = "Metric", values_from = "value") %>%
    bind_cols(nd) %>%
    mutate(Response = "Unweighted small",
    Management = ifelse(Management == "High", "High intensity", "Low intensity"),
    Sdep = 3.818393*Sdep - 4.769306)
pl_1 <- plot_dat %>%
  ggplot() +
  # geom_point(aes(x = Sdep, y = PH, colour = Management),
  #            alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_ribbon(data = f2,
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, 
                  group = Management),
              stat = "identity", 
              alpha = 1/4) +
  geom_smooth(data = filter(f, Response == "Unweighted small"),
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity", 
              alpha = 1/2, size = 1/2) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0))


# Plot for comparing Ndep effects on pH change
nd <- 
  tibble(Ndep = seq(from = -1.8, to = 3.15, length.out = 30) %>% 
           rep(., times = 2),
         fieldseason_rain = 0,
         Management = rep(c("Low","High"), each = 30),
         Year1_pH = 0,
         Sdep = 0,
         PH_SE = 0.3,
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted small")
)) %>%
  mutate(
    Ndep = 70.08984*Ndep + 138.8274
  )

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Year1_pH, Ndep,
         Sdep, fieldseason_rain) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))


plot_dat %>%
  ggplot() +
  geom_point(aes(x = Ndep, y = PH, colour = Management),
             alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity", 
              alpha = 1/2, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "N deposition", y = "pH change") +
  scale_x_continuous(expand = expansion(),limits = c(0,390)) 
ggsave("Ndep effect on pH with data measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/", 
       width = 16, height = 10, units = "cm")


# mean and prediction error - small unweighted only
f2 <- posterior_predict(mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
  as_tibble() %>%
  summarise(across(.fns = list(Estimate = mean,
                               Q2.5 = ~quantile(.x, probs = 0.025),
                               Q97.5 = ~quantile(.x, probs = 0.975)))) %>%
  pivot_longer(everything(), names_to = c("ID","Metric"),
               names_sep = "_") %>%
  pivot_wider(names_from = "Metric", values_from = "value") %>%
  bind_cols(nd) %>%
  mutate(Response = "Unweighted small",
         Management = ifelse(Management == "High", "High intensity", "Low intensity"),
         Ndep = 70.08984*Ndep + 138.8274)
pl_2 <- plot_dat %>%
  ggplot() +
  # geom_point(aes(x = Sdep, y = PH, colour = Management),
  #            alpha = 0.1) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_ribbon(data = f2,
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, 
                  group = Management),
              stat = "identity", 
              alpha = 1/4) +
  geom_smooth(data = filter(f, Response == "Unweighted small"),
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity", 
              alpha = 1/2, size = 1/2) +
  labs(x = "N deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0))

# Plot for comparing pH change on Ellenberg R change
nd <- 
  tibble(PH = seq(from = -3, to = 3.5, length.out = 40) %>% 
           rep(., times = 2),
         Management = rep(c("High","Low"), each = 40),
         Year1_pH = 0,
         N = 0,
         ELL_SE = 0.4,
         YRnm = 1,
         REP_ID = 1:80)

f <- do.call(rbind, list(
  fitted(mod_whw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(mod_whuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(mod_smw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(mod_smuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == "High", "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Year1_pH, N = NC_RATIO) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"))

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  scale_colour_manual(values = mang_cols,
                      aesthetics = c("colour","fill")) +
  facet_grid(~Response) +
  labs(x = "pH change", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("pH change effect on Ellenberg R measerror.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/", 
       width = 20, height = 10, units = "cm")

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = PH, y = value, colour = Management),
             alpha = 0.5) +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(~Response) +
  scale_colour_manual(values = mang_cols,
                      aesthetics = c("colour","fill")) +
  labs(x = "pH change", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("pH change effect on Ellenberg R with data measerror.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models", 
       width = 20, height = 10, units = "cm")
plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = PH, y = value, colour = Management),
             alpha = 0.5) +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_wrap(~Response) +
  scale_colour_manual(values = mang_cols,
                      aesthetics = c("colour","fill")) +
  labs(x = "pH change", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("pH change effect on Ellenberg R with data measerror.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models", 
       width = 15, height = 12, units = "cm")

# Parameter estimates
mcmc_plot(mod_whw, type = "areas_ridges") + 
  ggtitle(expression("Cover weighted Ellenberg R 200m"^2)) + 
  mcmc_plot(mod_whuw, type = "areas_ridges") + theme(axis.text.y = element_blank())+
  ggtitle(expression("Unweighted Ellenberg R 200m"^2)) +
  mcmc_plot(mod_smw, type = "areas_ridges") + 
  ggtitle(expression("Cover weighted Ellenberg R 4m"^2)) +
  mcmc_plot(mod_smuw, type = "areas_ridges") + theme(axis.text.y = element_blank()) +
  ggtitle(expression("Unweighted Ellenberg R 4m"^2)) 

ggsave("Parameter estimates bayesplot.png", 
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/",
       width = 24, height = 14, units = "cm", scale = 1.2)


mcmc_plot(mod_whw, type = "areas_ridges") +
  ggtitle(expression("Cover weighted Ellenberg R 200m"^2)) +
  mcmc_plot(mod_whuw, type = "areas_ridges") +
  ggtitle(expression("Unweighted Ellenberg R 200m"^2)) +
  mcmc_plot(mod_smw, type = "areas_ridges") +
  ggtitle(expression("Cover weighted Ellenberg R 4m"^2)) +
  mcmc_plot(mod_smuw, type = "areas_ridges") +
  ggtitle(expression("Unweighted Ellenberg R 4m"^2))  + plot_layout(ncol = 1)

ggsave("Parameter estimates bayesplot tall.png", 
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/",
       width = 10, height = 20, units = "cm", scale = 1.5)

plot_pars <- c(parnames(mod_whw)[1:21],"ar_PH","ar_Ell")
param_summaries <- do.call(rbind, list(
  posterior_summary(mod_whw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "WH_W"),
  posterior_summary(mod_whuw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "WH_UW"),
  posterior_summary(mod_smw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "SM_W"),
  posterior_summary(mod_smuw, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "SM_UW")
)) 
write.csv(param_summaries, 
          "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/Parameter_summary_table.csv",
          row.names = FALSE)


param_summaries <- param_summaries %>% 
  mutate(across(Estimate:Q97.5, round, 3)) %>%
  mutate(CI = paste0(Q2.5, " - ", Q97.5),
         Estimate = as.character(Estimate)) %>%
  select(-Est.Error, -Q2.5, -Q97.5) %>%
  pivot_longer(c(Estimate,CI)) %>%
  pivot_wider(names_from = c(Model, name), names_sep = "__") %>%
  mutate(across(ends_with("Estimate"),as.numeric))
writexl::write_xlsx(param_summaries, 
                    "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/Parameter_summary_nicetable.xlsx")


post <- do.call(rbind, list(
  posterior_samples(mod_whw, pars = plot_pars) %>% 
    mutate(Model = "WH_W"),
  posterior_samples(mod_whuw, pars = plot_pars) %>% 
    mutate(Model = "WH_UW"),
  posterior_samples(mod_smw, pars = plot_pars) %>% 
    mutate(Model = "SM_W"),
  posterior_samples(mod_smuw, pars = plot_pars) %>% 
    mutate(Model = "SM_UW")
)) 
post %>%
  transmute(Ell_pH_Low    = bsp_Ell_miPH + `bsp_Ell_miPH:ManagementLow`,
            Ell_pH_High = bsp_Ell_miPH,
            pH_Sdep_Low = b_PH_Sdep + `b_PH_ManagementLow:Sdep`,
            pH_Sdep_High = b_PH_Sdep,
            pH_Ndep_Low = b_PH_Ndep + `b_PH_ManagementLow:Ndep`,
            pH_Ndep_High = b_PH_Ndep,
            # pH_InitialpH_NA = b_PH_Year1_pH,
            # pH_Rainfall_NA = b_PH_fieldseason_rain,
            Model = Model) %>%
  tidyr::pivot_longer(contains("_"), "key", "value") %>%
  tidyr::separate(key, c("Response","Predictor","Management"),
                  sep = "_", remove = FALSE) %>%
  mutate(group = paste0(Model,key),
         Response = paste("Response:",Response),
         Predictor = paste("Predictor:",Predictor),
         Model = recode(Model,
                        "SM_UW" = "Small unweighted",
                        "SM_W" = "Small weighted",
                        "WH_W" = "Whole weighted",
                        "WH_UW" = "Whole unweighted")) %>%
  ggplot(aes(x = value, group = group, 
             color = Management, fill = Management)) +
  geom_density(alpha = 1/4) +
  facet_grid(Model~Response + Predictor) +
  scale_x_continuous("Parameter estimate", expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0,0)) +
  scale_colour_manual(values = mang_cols, aesthetics = c("colour","fill")) 
ggsave("Parameter estimates of pH and Ndep and Sdep by management.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models",
       width = 15, height = 15, units = "cm")

# NS as Ellenberg predictors
mcmc_plot(mod_whw_ns, type = "areas_ridges") + 
  ggtitle(expression("Cover weighted Ellenberg R 200m"^2)) + 
  mcmc_plot(mod_whuw_ns, type = "areas_ridges") + theme(axis.text.y = element_blank())+
  ggtitle(expression("Unweighted Ellenberg R 200m"^2)) +
  mcmc_plot(mod_smw_ns, type = "areas_ridges") + 
  ggtitle(expression("Cover weighted Ellenberg R 4m"^2)) +
  mcmc_plot(mod_smuw_ns, type = "areas_ridges") + theme(axis.text.y = element_blank()) +
  ggtitle(expression("Unweighted Ellenberg R 4m"^2)) 

ggsave("Parameter estimates Ell_NS models bayesplot.png", 
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/",
       width = 24, height = 14, units = "cm", scale = 1.3)


plot_pars <- c(parnames(mod_whw_ns)[1:25],"ar_PH","ar_Ell")
param_summaries <- do.call(rbind, list(
  posterior_summary(mod_whw_ns, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "WH_W"),
  posterior_summary(mod_whuw_ns, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "WH_UW"),
  posterior_summary(mod_smw_ns, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "SM_W"),
  posterior_summary(mod_smuw_ns, pars = plot_pars) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Parameter") %>%
    mutate(Model = "SM_UW")
)) %>% 
  mutate(across(Estimate:Q97.5, round, 3)) %>%
  mutate(CI = paste0(Q2.5, " - ", Q97.5),
         Estimate = as.character(Estimate)) %>%
  select(-Est.Error, -Q2.5, -Q97.5) %>%
  pivot_longer(c(Estimate,CI)) %>%
  pivot_wider(names_from = c(Model, name), names_sep = "__") %>%
  mutate(across(ends_with("Estimate"),as.numeric))
writexl::write_xlsx(param_summaries, 
                    "Outputs/Models/Difference/Multivariate_measerrorXYU/Linear_models/Parameter_summary_nsmodel_nicetable.xlsx")

