# Ellenberg R and pH model
library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

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
## ** Sdep ####
# model each Ellenberg R change as a function of Sdep
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W, Ell_SE = ELL_WH_W_SE_NORM, 
         Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ell_SE, Sdep, Management) %>%
  na.omit()

# prior checks
get_prior(Ell | mi(Ell_SE) ~ Sdep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
          data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(normal(0.4,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))

# prior predictive check model
prior_mod <- brm(Ell |mi(Ell_SE) ~ Sdep*Management + 
                   (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                 sample_prior = "only",
                 data = mod_data, prior = mod_pr, cores = 4, iter = 5000)
summary(prior_mod)

# ellenberg R model - weighted Ellenberg R for whole plot
ell_Sdep_diff <- brm(Ell | mi(Ell_SE) ~ Sdep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                     data = mod_data, prior = mod_pr, 
                     cores = 4, iter = 5000, 
                     file = "Outputs/Models/Difference/Univariate_measerror/EllR_WHW_Sdep_HAB")
summary(ell_Sdep_diff)
plot(ell_Sdep_diff, ask = FALSE)
pp_check(ell_Sdep_diff)

# unweighted ellenberg r whole plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW, Ell_SE = ELL_WH_UW_SE_NORM,
         Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ell_SE, Sdep, Management) %>%
  na.omit() 
ell_Sdep_diffuw <- update(ell_Sdep_diff, newdata = mod_data, cores=4, iter = 5000,
                          file = "Outputs/Models/Difference/Univariate_measerror/EllR_WHUW_Sdep_HAB")
summary(ell_Sdep_diffuw)
plot(ell_Sdep_diffuw, ask = FALSE)
pp_check(ell_Sdep_diffuw)


# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W, Ell_SE = ELL_SM_W_SE_NORM,
         Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ell_SE, Sdep, Management) %>%
  na.omit() 
ell_Sdep_diffsmw <- update(ell_Sdep_diff, newdata = mod_data, cores=4, iter = 5000,
                           control = list(adapt_delta = 0.99),
                           file = "Outputs/Models/Difference/Univariate_measerror/EllR_SMW_Sdep_HAB")
summary(ell_Sdep_diffsmw)
plot(ell_Sdep_diffsmw, ask = FALSE)
pp_check(ell_Sdep_diffsmw)

# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW, Ell_SE = ELL_SM_UW_SE_NORM,
         Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ell_SE, Sdep, Management) %>%
  na.omit() 
ell_Sdep_diffsmuw <- update(ell_Sdep_diff, newdata = mod_data, cores=4, iter = 5000,
                            control = list(adapt_delta = 0.99),
                            file = "Outputs/Models/Difference/Univariate_measerror/EllR_SMUW_Sdep_HAB")
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

# ** Ndep ####
# model each Ellenberg R change as a function of Ndep
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W, Ell_SE = ELL_WH_W_SE_NORM,
         Management = ifelse(Management == "High",1,0),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ell_SE, Ndep, Management) %>%
  na.omit()

# prior checks
get_prior(Ell | mi(Ell_SE) ~ Ndep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
          family = "student", data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(normal(0.4,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))

# prior predictive check model
prior_mod <- brm(Ell | mi(Ell_SE) ~ Ndep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                 family = "student", sample_prior = "only",
                 data = mod_data, prior = mod_pr, cores = 4, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50)
# underestimates peak at 0 on average but gets close

# ellenberg R model - weighted Ellenberg R for whole plot
ell_Ndep_diff <- brm(Ell | mi(Ell_SE) ~ Ndep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                     data = mod_data, prior = mod_pr, family = "student",
                     cores = 4, iter = 5000, 
                     file = "Outputs/Models/Difference/Univariate_measerror/EllR_WHW_Ndep_HAB")
summary(ell_Ndep_diff)
plot(ell_Ndep_diff, ask = FALSE)
pp_check(ell_Ndep_diff)

# unweighted ellenberg r whole plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW, Ell_SE = ELL_WH_UW_SE_NORM,
         Management = ifelse(Management == "High",1,0),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ell_SE, Ndep, Management) %>%
  na.omit() 
ell_Ndep_diffuw <- update(ell_Ndep_diff, newdata = mod_data, cores=4, iter = 5000,
                          file = "Outputs/Models/Difference/Univariate_measerror/EllR_WHUW_Ndep_HAB")
summary(ell_Ndep_diffuw)
plot(ell_Ndep_diffuw, ask = FALSE)
pp_check(ell_Ndep_diffuw)


# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W, Ell_SE = ELL_SM_W_SE_NORM,
         Management = ifelse(Management == "High",1,0),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ell_SE, Ndep, Management) %>%
  na.omit() 
ell_Ndep_diffsmw <- update(ell_Ndep_diff, newdata = mod_data, cores = 4, iter = 5000,
                           control = list(adapt_delta = 0.99),
                           file = "Outputs/Models/Difference/Univariate_measerror/EllR_SMW_Ndep_HAB")
summary(ell_Ndep_diffsmw)
plot(ell_Ndep_diffsmw, ask = FALSE)
pp_check(ell_Ndep_diffsmw)

# unweighted ellenberg r small plot
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW, Ell_SE = ELL_SM_UW_SE_NORM,
         Management = ifelse(Management == "High",1,0),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Ell_SE, Ndep, Management) %>%
  na.omit() 
ell_Ndep_diffsmuw <- update(ell_Ndep_diff, newdata = mod_data, cores = 4, iter = 5000,
                            control = list(adapt_delta = 0.99),
                            file = "Outputs/Models/Difference/Univariate_measerror/EllR_SMUW_Ndep_HAB")
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


# *** pH ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Management = ifelse(Management == "High",1,0),
         PH_SE = PH_DIW_SE_NORM,
         Sdep = as.numeric(scale(Sdep)),
         Ndep = as.numeric(scale(Ndep))) %>%
  select(REP_ID, SQUARE, YRnm, PH, PH_SE, Sdep, Ndep, Management) %>%
  na.omit() 
ph_Sdep_diff <- brm(PH | mi(PH_SE) ~ Sdep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                    data = mod_data, prior = mod_pr,
                    cores = 4, iter = 5000, 
                    file = "Outputs/Models/Difference/Univariate_measerror/PH_Sdep_HAB")
summary(ph_Sdep_diff)
plot(ph_Sdep_diff, ask = FALSE)
pp_check(ph_Sdep_diff)

ph_Ndep_diff <- brm(PH | mi(PH_SE) ~ Ndep*Management + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                    data = mod_data, prior = mod_pr,
                    cores = 4, iter = 5000, 
                    file = "Outputs/Models/Difference/Univariate_measerror/PH_Ndep_HAB")
summary(ph_Ndep_diff)
plot(ph_Ndep_diff, ask = FALSE)
pp_check(ph_Ndep_diff)


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

util <- new.env()
source('stan_utility.R', local=util)
c_dark_trans <- "#80808080"
c_yellow_trans <- "#FFFF0080"

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
summary(full_mod_hab)
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
get_prior(bf(Ell  ~ Management*PH + s(Year1_pH, N, Management) + 
               (1|SQUARE) +
               ar(time = YRnm, gr = REP_ID)) +
            bf(PH  ~ Management*Sdep + fieldseason_rain + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID)) + 
            bf(N ~ Ndep + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID)) +
            set_rescor(FALSE), data = mod_data)


mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.5), class = "b", resp = "N"),
            prior(student_t(3, 0, 2.5), class = "sds", resp = "Ell"),
            prior(student_t(3, 0, 2.5), class = "sds", resp = "PH"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(student_t(3,0,1), class = "Intercept", resp = "N"),
            prior(normal(0,0.2), class = "ar", resp = "Ell"),
            prior(normal(0,0.2), class = "ar", resp = "PH"),
            prior(normal(0,0.2), class = "ar", resp = "N"),
            prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sd", resp = "N"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "N"))


# prior simulation
prior_mod <- brm(bf(Ell  ~ Management*PH + s(Year1_pH, N, Management) + 
                      (1|SQUARE) +
                      ar(time = YRnm, gr = REP_ID)) +
                   bf(PH  ~ Management*Sdep + fieldseason_rain  + s(Year1_pH, N, Management) +
                        (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID)) + 
                   bf(N ~ Ndep + (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID))  +
                   set_rescor(FALSE), data = mod_data, prior = mod_pr,
                 sample_prior = "only", save_pars = save_pars(all = TRUE, latent = TRUE), 
                 cores = 4, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50, resp = "Ell")
pp_check(prior_mod, nsamples = 50, resp = "PH")
pp_check(prior_mod, nsamples = 50, resp = "N") 
plot(conditional_effects(prior_mod))

# run model - full weighted Ell R
full_mod_whw <- brm(bf(Ell  ~ Management*PH + s(Year1_pH, N) + 
                         (1|SQUARE) +
                         ar(time = YRnm, gr = REP_ID)) +
                      bf(PH  ~ Management*Sdep + fieldseason_rain  + 
                           s(Year1_pH, N) +
                           (1|SQUARE) +
                           ar(time = YRnm, gr = REP_ID)) + 
                      bf(N ~ Ndep + (1|SQUARE) +
                           ar(time = YRnm, gr = REP_ID)) +
                      set_rescor(FALSE), data = mod_data, prior = mod_pr,
                    file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHW_spl_PH_spl_N_HAB",
                    cores = 4, iter = 4000, control = list(adapt_delta = 0.95))
summary(full_mod_whw)
plot(full_mod_whw, ask = FALSE)
pp_check(full_mod_whw, nsamples = 50, resp = "Ell")
pp_check(full_mod_whw, nsamples = 50, resp = "PH")
pp_check(full_mod_whw, nsamples = 50, resp = "N")

# plot(conditional_effects(full_mod_whw))

# No pH ~ N spline
mod_pr_red1 <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
                 prior(normal(0,0.5), class = "b", resp = "PH"),
                 prior(normal(0,0.5), class = "b", resp = "N"),
                 prior(student_t(3, 0, 2.5), class = "sds", resp = "Ell"),
                 prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
                 prior(normal(0,0.25), class = "Intercept", resp = "PH"),
                 prior(student_t(3,0,1), class = "Intercept", resp = "N"),
                 prior(normal(0,0.2), class = "ar", resp = "Ell"),
                 prior(normal(0,0.2), class = "ar", resp = "PH"),
                 prior(normal(0,0.2), class = "ar", resp = "N"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "N"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "N"))

red_mod_1_whw <- brm(bf(Ell  ~ Management*PH + s(Year1_pH, N) + 
                          (1|SQUARE) +
                          ar(time = YRnm, gr = REP_ID),
                        family = "student") +
                       bf(PH  ~ Management*Sdep + fieldseason_rain +
                            (1|SQUARE) +
                            ar(time = YRnm, gr = REP_ID), family = "student") + 
                       bf(N ~ Ndep + (1|SQUARE) +
                            ar(time = YRnm, gr = REP_ID)) +
                       set_rescor(FALSE), data = mod_data, prior = mod_pr_red1,
                     file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHW_spl_PH_N_HAB",
                     cores = 4, iter = 4000, control = list(adapt_delta = 0.99))
summary(red_mod_1_whw)
plot(red_mod_1_whw, ask = FALSE)
pp_check(red_mod_1_whw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_whw, nsamples = 50, resp = "PH")
pp_check(red_mod_1_whw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_1_whw))

# No Ell ~ N spline
mod_pr_red2 <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
                 prior(normal(0,0.5), class = "b", resp = "PH"),
                 prior(normal(0,0.5), class = "b", resp = "N"),
                 prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
                 prior(normal(0,0.25), class = "Intercept", resp = "PH"),
                 prior(student_t(3,0,1), class = "Intercept", resp = "N"),
                 prior(normal(0,0.2), class = "ar", resp = "Ell"),
                 prior(normal(0,0.2), class = "ar", resp = "PH"),
                 prior(normal(0,0.2), class = "ar", resp = "N"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "N"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "N"))
red_mod_2_whw <- brm(bf(Ell  ~ Management*PH + 
                          (1|SQUARE) +
                          ar(time = YRnm, gr = REP_ID)) +
                       bf(PH  ~ Management*Sdep + fieldseason_rain +
                            (1|SQUARE) +
                            ar(time = YRnm, gr = REP_ID)) + 
                       bf(N ~ Ndep + (1|SQUARE) +
                            ar(time = YRnm, gr = REP_ID)) +
                       set_rescor(FALSE), data = mod_data, prior = mod_pr_red2,
                     file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHW_PH_N_HAB",
                     cores = 4, iter = 4000, control = list(adapt_delta = 0.95))
summary(red_mod_2_whw)
plot(red_mod_2_whw, ask = FALSE)
pp_check(red_mod_2_whw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_whw, nsamples = 50, resp = "PH")
pp_check(red_mod_2_whw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_2_whw))

loo_compare(full_mod_whw, red_mod_1_whw, red_mod_2_whw, criterion = "kfold")
# elpd_diff se_diff
# red_mod_1_whw    0.0       0.0 
# full_mod_whw   -24.3      33.5 
# red_mod_2_whw -163.3      31.5  


# ~~~ 200m2 unweighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW) %>%
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
# run model - full unweighted Ell R
full_mod_whuw <- update(full_mod_whw, newdata = mod_data,
                        control = list(adapt_delta = 0.95),
                        cores = 4, iter = 4000,
                        file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHUW_spl_PH_spl_N_HAB")

summary(full_mod_whuw)
plot(full_mod_whuw, ask = FALSE)
pp_check(full_mod_whuw, nsamples = 50, resp = "Ell")
pp_check(full_mod_whuw, nsamples = 50, resp = "PH")
pp_check(full_mod_whuw, nsamples = 50, resp = "N")

# plot(conditional_effects(full_mod_whuw))


# No pH ~ N spline
red_mod_1_whuw <- update(red_mod_1_whw, newdata = mod_data,
                         cores = 4, iter = 4000,
                         control = list(adapt_delta = 0.99),
                         file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHUW_spl_PH_N_HAB")

summary(red_mod_1_whuw)
plot(red_mod_1_whuw, ask = FALSE)
pp_check(red_mod_1_whuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_whuw, nsamples = 50, resp = "PH")
pp_check(red_mod_1_whuw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_1_whuw))

# No Ell ~ N spline
red_mod_2_whuw <- update(red_mod_2_whw, newdata = mod_data,
                         cores = 4, iter = 4000,
                         control = list(adapt_delta = 0.95),
                         file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHUW_PH_N_HAB")

summary(red_mod_2_whuw)
plot(red_mod_2_whuw, ask = FALSE)
pp_check(red_mod_2_whuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_whuw, nsamples = 50, resp = "PH")
pp_check(red_mod_2_whuw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_2_whuw))
loo_compare(full_mod_whuw, red_mod_1_whuw, red_mod_2_whuw,
            criterion = "kfold")
# elpd_diff se_diff
# red_mod_1_whuw    0.0       0.0 
# full_mod_whuw  -127.3      44.2 
# red_mod_2_whuw -187.3      44.4

# ~~~ 4m2 weighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W) %>%
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
# run model - small weighted Ell R
full_mod_smw <- update(full_mod_whw, newdata = mod_data,
                       cores = 4, iter = 4000,
                       control = list(adapt_delta = 0.95),
                       file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMW_spl_PH_spl_N_HAB")

summary(full_mod_smw)
plot(full_mod_smw, ask = FALSE)
pp_check(full_mod_smw, nsamples = 50, resp = "Ell")
pp_check(full_mod_smw, nsamples = 50, resp = "PH")
pp_check(full_mod_smw, nsamples = 50, resp = "N")

# plot(conditional_effects(full_mod_smw))


# No pH ~ N spline
red_mod_1_smw <- update(red_mod_1_whw, newdata = mod_data,
                        cores = 4, iter = 4000,
                        control = list(adapt_delta = 0.99),
                        file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMW_spl_PH_N_HAB")

summary(red_mod_1_smw)
plot(red_mod_1_smw, ask = FALSE)
pp_check(red_mod_1_smw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_smw, nsamples = 50, resp = "PH")
pp_check(red_mod_1_smw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_1_smw))

# No Ell ~ N spline
red_mod_2_smw <- update(red_mod_2_whw, newdata = mod_data,
                        cores = 4, iter = 4000,
                        control = list(adapt_delta = 0.95),
                        file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMW_PH_N_HAB")

summary(red_mod_2_smw)
plot(red_mod_2_smw, ask = FALSE)
pp_check(red_mod_2_smw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_smw, nsamples = 50, resp = "PH")
pp_check(red_mod_2_smw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_2_smw))
loo_compare(full_mod_smw, red_mod_1_smw, red_mod_2_smw,
            criterion = "kfold")
# elpd_diff se_diff
# full_mod_smw     0.0       0.0 
# red_mod_1_smw   -5.9      33.3 
# red_mod_2_smw -151.0      20.4

# ~~~ 4m2 unweighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW) %>%
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


# Multivariate with N as a 2D spline
# Habitat interaction
# run model - small unweighted Ell R
full_mod_smuw <- update(full_mod_whw, newdata = mod_data,
                        cores = 4, iter = 4000,
                        control = list(adapt_delta = 0.95),
                        file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMUW_spl_PH_spl_N_HAB")

summary(full_mod_smuw)
plot(full_mod_smuw, ask = FALSE)
pp_check(full_mod_smuw, nsamples = 50, resp = "Ell")
pp_check(full_mod_smuw, nsamples = 50, resp = "PH")
pp_check(full_mod_smuw, nsamples = 50, resp = "N")

# mod_pr_norm <- c(prior(normal(0,0.5), class = "b"),
#                  prior(student_t(3, 0, 2.5), class = "sds"),
#                  prior(normal(0,0.25), class = "Intercept"),
#                  prior(normal(0,0.2), class = "ar"),
#                  prior(student_t(3,0,0.5), class = "sd"),
#                  prior(student_t(3,0,0.5), class = "sigma"))
# 
# test_mod_norm <- brm(Ell  ~ Management*PH + s(Year1_pH, N, Management) + 
#                       (1|SQUARE) + ar(time = YRnm, gr = REP_ID), 
#                     data = mod_data, 
#                     prior = mod_pr_norm,
#                     save_pars = save_pars(all = TRUE), 
#                     cores = 4, iter = 5000)

# plot(conditional_effects(full_mod_smuw))


# No pH ~ N spline
red_mod_1_smuw <- update(red_mod_1_whw, newdata = mod_data,
                         cores = 4, iter = 4000,
                         control = list(adapt_delta = 0.99),
                         file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMUW_spl_PH_N_HAB")

summary(red_mod_1_smuw)
plot(red_mod_1_smuw, ask = FALSE)
pp_check(red_mod_1_smuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_smuw, nsamples = 50, resp = "PH")
pp_check(red_mod_1_smuw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_1_smuw))

# No Ell ~ N spline
red_mod_2_smuw <- update(red_mod_2_whw, newdata = mod_data,
                         cores = 4, iter = 4000,
                         control = list(adapt_delta = 0.95),
                         file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMUW_PH_N_HAB")

summary(red_mod_2_smuw)
plot(red_mod_2_smuw, ask = FALSE)
pp_check(red_mod_2_smuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_smuw, nsamples = 50, resp = "PH")
pp_check(red_mod_2_smuw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_2_smuw))

loo_compare(full_mod_smuw, red_mod_1_smuw, red_mod_2_smuw,
            criterion = "kfold")
# elpd_diff se_diff
# red_mod_1_smuw    0.0       0.0 
# full_mod_smuw   -44.1      37.7 
# red_mod_2_smuw -200.8      33.1 



# Summary plots ####
# Plot for comparing initial pH effects on pH change
nd <- 
  tibble(Year1_pH = seq(from = -2.25, to = 3.6, length.out = 30) %>% 
           rep(., times = 6),
         N = rep(c(-1,0,1), each = 30) %>%
           rep(., times = 2),
         Management = rep(0:1, each = 90),
         Sdep = 0,
         fieldseason_rain = 0,
         YRnm = 1,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
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
  mutate(N = as.character(N)) %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)) %>%
                filter(Management == "Low intensity"),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Initial pH", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("Initial pH and N effect on pH under low intensity.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")

plot_dat %>%
  mutate(N = cut(scale(N), c(-10,-0.5,0.5,10), labels = c("-1","0","1")),
         Year1_pH = scale(Year1_pH)) %>%
  filter(!is.na(N)) %>%
  ggplot() +
  geom_point(aes(x = Year1_pH, y = PH, colour = N),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "Initial pH", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Initial pH and N effect on pH with data facet by management.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")


nd <- 
  tibble(Year1_pH = seq(from = -2.25, to = 3.6, length.out = 30) %>% 
           rep(., times = 6),
         N = rep(c(-1,0,1), each = 30) %>%
           rep(., times = 2),
         Management = rep(0.5, each = 180),
         Sdep = 0,
         fieldseason_rain = 0,
         YRnm = 1,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat %>%
  mutate(N = cut(scale(N), c(-10,-0.5,0.5,10), labels = c("-1","0","1")),
         Year1_pH = scale(Year1_pH)) %>%
  filter(!is.na(N)) %>%
  ggplot() +
  geom_point(aes(x = Year1_pH, y = PH, colour = N),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)) %>%
                filter(Management == "Low intensity"),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Initial pH", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Initial pH and N effect on pH with data predict middle management.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 18, height = 12, units = "cm")

# Plot for comparing initial pH effects on pH change
nd <- 
  tibble(Sdep = seq(from = -5.5, to = 1.75, length.out = 30) %>% 
           rep(., times = 6),
         fieldseason_rain = rep(c(-1,0,1), each = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 90),
         Year1_pH = 0,
         N = 0,
         YRnm = 1,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

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
         Management = ifelse(Management == "High", "High intensity", "Low intensity"),
         Sdep =scale(Sdep),
         fieldseason_rain = scale(fieldseason_rain))

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, fieldseason_rain = as.character(fieldseason_rain)),
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = fieldseason_rain, color = fieldseason_rain,
                  group = fieldseason_rain),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Sdep and rain effect on pH.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")

plot_dat %>%
  mutate(fieldseason_rain = cut(scale(fieldseason_rain), c(-10,-0.5,0.5,10), labels = c("-1","0","1"))) %>%
  filter(!is.na(fieldseason_rain)) %>%
  ggplot() +
  geom_point(aes(x = Sdep, y = PH, colour = fieldseason_rain),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, fieldseason_rain = as.character(fieldseason_rain)),
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = fieldseason_rain, color = fieldseason_rain,
                  group = fieldseason_rain),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Sdep and rain effect on pH with data.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")

# Plot for comparing initial pH effects on Ellenberg R change
nd <- 
  tibble(Year1_pH = seq(from = -2.25, to = 3.6, length.out = 30) %>% 
           rep(., times = 6),
         N = rep(c(-1,0,1), each = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 90),
         Sdep = 0,
         PH = 0,
         fieldseason_rain = 0,
         YRnm = 1,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
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
  mutate(N = as.character(N)) %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "Initial pH", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("Initial pH and N effect on Ellenberg R.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")

plot_dat %>%
  mutate(N = cut(scale(N), c(-10,-0.5,0.5,10), labels = c("-1","0","1")),
         Year1_pH = scale(Year1_pH)) %>%
  filter(!is.na(N)) %>%
  ggplot() +
  geom_point(aes(x = Year1_pH, y = value, colour = N),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "Initial pH", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Initial pH and N effect on Ellenberg R with data.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")

# Plot for comparing pH change on Ellenberg R change
nd <- 
  tibble(PH = seq(from = -3, to = 3.5, length.out = 40) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 40),
         Year1_pH = 0,
         N = 0,
         YRnm = 1,
         REP_ID = 1:80)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
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
ggsave("pH change effect on Ellenberg R.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 20, height = 10, units = "cm")

plot_dat %>%
  filter(!is.na(N)) %>%
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
ggsave("pH change effect on Ellenberg R with data.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 20, height = 10, units = "cm")

# ~~ no splines ####
# Plot for comparing pH change on Ellenberg R change
nd <- 
  tibble(PH = seq(from = -3, to = 3.5, length.out = 40) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 40),
         YRnm = 1,
         REP_ID = 1:80)

f <- do.call(rbind, list(
  fitted(red_mod_2_whw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(red_mod_2_whuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(red_mod_2_smw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(red_mod_2_smuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
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
ggsave("pH change effect on Ellenberg R no spline multivariate model.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 20, height = 10, units = "cm")

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = PH, y = value, colour = Management),
             alpha = 1/3) +
  geom_smooth(data = f,
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/2, size = 1/2) +
  facet_grid(~Response) +
  labs(x = "pH change", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("pH change effect on Ellenberg R with data no spline model.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 20, height = 10, units = "cm")

# Plot for comparing initial pH effects on pH change
nd <- 
  tibble(Sdep = seq(from = -5.5, to = 1.75, length.out = 30) %>% 
           rep(., times = 6),
         fieldseason_rain = rep(c(-1,0,1), each = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 90),
         YRnm = 1,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(red_mod_2_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(red_mod_2_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(red_mod_2_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(red_mod_2_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Sdep, fieldseason_rain) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"),
         Sdep =scale(Sdep),
         fieldseason_rain = scale(fieldseason_rain))

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, fieldseason_rain = as.character(fieldseason_rain)),
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = fieldseason_rain, color = fieldseason_rain,
                  group = fieldseason_rain),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Sdep and rain effect on pH no spline model.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")

plot_dat %>%
  mutate(fieldseason_rain = cut(scale(fieldseason_rain), c(-10,-0.5,0.5,10), labels = c("-1","0","1"))) %>%
  filter(!is.na(fieldseason_rain)) %>%
  ggplot() +
  geom_point(aes(x = Sdep, y = PH, colour = fieldseason_rain),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, fieldseason_rain = as.character(fieldseason_rain)),
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = fieldseason_rain, color = fieldseason_rain,
                  group = fieldseason_rain),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Sdep and rain effect on pH no spline model with data.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")

# Plot for comparing impact of Ndep upon N:C
nd <- 
  tibble(Ndep = seq(from = -1.8, to = 3.2, length.out = 60),
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "N") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "N") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "N") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "N") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  filter(!is.na(Management)) %>%
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         N = NC_RATIO, Ndep, REP_ID, Time_period) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Ndep =scale(Ndep),
         N = scale(N))

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity", 
              alpha = 1/4, size = 1/2,
              fill = "#CC79A7", colour = "#CC79A7") +
  facet_wrap(~Response) +
  labs(x = "N deposition", y = "N:C ratio") +
  scale_x_continuous(expand = c(0,0))
ggsave("Ndep effect on NC full model.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")

plot_dat %>%
  ggplot() +
  geom_point(aes(x = Ndep, y = N),
             alpha = 0.5, colour = "#CC79A7") +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity", 
              alpha = 1/3, size = 1/2,
              fill = "black", colour = "black") +
  facet_wrap(~Response) +
  labs(x = "N deposition", y = "N:C ratio") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("Ndep effect on NC full model with data.png",
       path = "Outputs/Models/Difference/Multivariate_nomeaserror", 
       width = 16, height = 18, units = "cm")


# ** Measurement error ####
# Multivariate with Ndep as a 2D spline
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
  mutate(Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)),
         N = as.numeric(scale(N)),
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()


# Habitat interaction
get_prior(bf(Ell | mi(ELL_SE) ~ Management*mi(PH) + s(Year1_pH, Ndep, Management) + 
               (1|p|SQUARE) +
               ar(time = YRnm, gr = REP_ID)) +
            bf(PH | mi(PH_SE)  ~ Management*Sdep + fieldseason_rain + s(Year1_pH, Ndep, Management) +
                 (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID)) + 
            set_rescor(FALSE), data = mod_data)


mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(student_t(3, 0, 2.5), class = "sds", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(normal(0,0.2), class = "ar", resp = "Ell"),
            prior(normal(0,0.2), class = "ar", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
            prior(lkj(5), class = "cor"))

mod_data_low <- filter(mod_data, Management == 0)
mod_pr_test <- c(prior(normal(0,0.5), class = "b"),
                 prior(student_t(3, 0, 2.5), class = "sds"),
                 prior(normal(0,0.25), class = "Intercept"),
                 prior(normal(0,0.2), class = "ar"),
                 prior(student_t(3,0,0.5), class = "sd"),
                 prior(student_t(3,0,0.5), class = "sigma"))
ph_test_mod <- brm(PH|mi(PH_SE) ~ fieldseason_rain + s(Year1_pH, Ndep, Sdep) + 
                     ar(time = YRnm, gr = REP_ID) + (1|SQUARE), 
                   data = mod_data_low, prior = mod_pr_test, 
                   cores = 4, iter = 4000, 
                   save_pars = save_pars(all = TRUE, latent = TRUE),
                   control = list(adapt_delta = 0.95))

ph_test_mod2 <- brm(PH|mi(PH_SE) ~ fieldseason_rain + Sdep + s(Year1_pH, Ndep) + 
                      ar(time = YRnm, gr = REP_ID) + (1|SQUARE), 
                    data = mod_data_low, prior = mod_pr_test, 
                    cores = 4, iter = 4000, 
                    save_pars = save_pars(all = TRUE, latent = TRUE))
ph_test_mod3 <- brm(PH|mi(PH_SE) ~ fieldseason_rain + Sdep + s(Year1_pH, N) + 
                      ar(time = YRnm, gr = REP_ID) + (1|SQUARE), 
                    data = mod_data_low, prior = mod_pr_test, 
                    cores = 4, iter = 4000, 
                    save_pars = save_pars(all = TRUE, latent = TRUE))


ell_test_mod2 <- brm(Ell|mi(ELL_SE) ~ fieldseason_rain + Sdep + s(Year1_pH, Ndep) + 
                       ar(time = YRnm, gr = REP_ID) + (1|SQUARE), 
                     data = mod_data_low, prior = mod_pr_test, 
                     cores = 4, iter = 4000, 
                     save_pars = save_pars(all = TRUE, latent = TRUE),
                     control = list(adapt_delta = 0.95))
ell_test_mod3 <- brm(Ell|mi(ELL_SE) ~ fieldseason_rain + Sdep + s(Year1_pH, N) + 
                       ar(time = YRnm, gr = REP_ID) + (1|SQUARE), 
                     data = mod_data_low, prior = mod_pr_test, 
                     cores = 4, iter = 4000, 
                     save_pars = save_pars(all = TRUE, latent = TRUE))


# prior simulation
prior_mod <- brm(bf(Ell | mi(ELL_SE) ~ Management*mi(PH) + s(Year1_pH, Ndep, Management) + 
                      (1|p|SQUARE) +
                      ar(time = YRnm, gr = REP_ID)) +
                   bf(PH | mi(PH_SE)  ~ Management*Sdep + fieldseason_rain + (1|p|SQUARE) +
                        ar(time = YRnm, gr = REP_ID)) + 
                   set_rescor(FALSE), data = mod_data, prior = mod_pr,
                 sample_prior = "only", save_pars = save_pars(all = TRUE, latent = TRUE), 
                 cores = 4, iter = 5000)
summary(prior_mod)
plot(prior_mod)
pp_check(prior_mod, nsamples = 50, resp = "Ell")
pp_check(prior_mod, nsamples = 50, resp = "PH")
pp_check(prior_mod, nsamples = 50, resp = "N") 
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
                                      Improved*0.5*pH_diffYear2 + 
                                      (pH>5.5)*-0.5*Ndep,EllSE_year1)) %>%
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
                                      Improved*0.5*pH_diffYear3 + 
                                      (pH_year2>5.5)*-0.5*Ndep,EllSE_year2)) %>%
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


sim_mod <- brm(bf(Ell | mi(ELL_SE) ~ Improved*mi(PH) + s(Year1_pH, Ndep, Improved) + 
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
  mutate(Management = ifelse(Management == "High",1,0),
         Sdep = as.numeric(scale(Sdep)), 
         Ndep = as.numeric(scale(Ndep)), 
         Year1_pH = as.numeric(scale(Year1_pH, scale = FALSE)),
         fieldseason_rain = as.numeric(scale(fieldseason_rain))) %>%
  na.omit()


# run model - full weighted Ell R
full_mod_whw_nocor <- brm(bf(Ell | mi(ELL_SE)  ~ Management*mi(PH) + s(Year1_pH, Ndep) + 
                               (1|SQUARE) +
                               ar(time = YRnm, gr = REP_ID)) +
                            bf(PH | mi(PH_SE)  ~ Management*Sdep + fieldseason_rain  + 
                                 s(Year1_pH, Ndep) +
                                 (1|SQUARE) +
                                 ar(time = YRnm, gr = REP_ID)) + 
                            set_rescor(FALSE), data = mod_data, prior = mod_pr[1:11,],
                          save_pars = save_pars(all = TRUE, latent = TRUE), 
                          # file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_WHW_spl_PH_spl_HAB",
                          cores = 4, iter = 4000)
summary(full_mod_whw)
plot(full_mod_whw, ask = FALSE)
pp_check(full_mod_whw, nsamples = 50, resp = "Ell")
pp_check(full_mod_whw, nsamples = 50, resp = "PH")

# plot(conditional_effects(full_mod_whw))
# check residuals against habitat
full_mod_whw_resid <- resid(full_mod_whw)
dimnames(full_mod_whw_resid)
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
                  Ell_resid = full_mod_whw_resid[,1,"Ell"],
                  PH_resid = full_mod_whw_resid[,1,"PH"]) %>%
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
ggsave("Residuals by habitat whole weighted full model.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")

# No pH ~ N spline
mod_pr_red1 <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
                 prior(normal(0,0.5), class = "b", resp = "PH"),
                 prior(student_t(3, 0, 2.5), class = "sds", resp = "Ell"),
                 prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
                 prior(normal(0,0.25), class = "Intercept", resp = "PH"),
                 prior(normal(0,0.2), class = "ar", resp = "Ell"),
                 prior(normal(0,0.2), class = "ar", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
                 prior(lkj(2), class = "cor"))

red_mod_1_whw <- brm(bf(Ell | mi(ELL_SE) ~ Management*mi(PH) + s(Year1_pH, Ndep) + 
                          (1|p|SQUARE) +
                          ar(time = YRnm, gr = REP_ID)) +
                       bf(PH | mi(PH_SE)  ~ Management*Sdep + fieldseason_rain +
                            (1|p|SQUARE) +
                            ar(time = YRnm, gr = REP_ID)) + 
                       set_rescor(FALSE), data = mod_data, prior = mod_pr_red1,
                     save_pars = save_pars(all = TRUE, latent = TRUE), 
                     file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_WHW_spl_PH_HAB",
                     cores = 4, iter = 4000, control = list(adapt_delta = 0.99))
summary(red_mod_1_whw)
plot(red_mod_1_whw, ask = FALSE)
pp_check(red_mod_1_whw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_whw, nsamples = 50, resp = "PH")

# plot(conditional_effects(red_mod_1_whw))
# check residuals against habitat
mod_resid <- resid(red_mod_1_whw)
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
ggsave("Residuals by habitat whole weighted reduced model 1.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")


# No Ell ~ N spline
mod_pr_red2 <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
                 prior(normal(0,0.5), class = "b", resp = "PH"),
                 prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
                 prior(normal(0,0.25), class = "Intercept", resp = "PH"),
                 prior(normal(0,0.2), class = "ar", resp = "Ell"),
                 prior(normal(0,0.2), class = "ar", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
                 prior(lkj(2), class = "cor"))
red_mod_2_whw <- brm(bf(Ell | mi(ELL_SE) ~ Management*mi(PH) + 
                          (1|p|SQUARE) +
                          ar(time = YRnm, gr = REP_ID)) +
                       bf(PH | mi(PH_SE)  ~ Management*Sdep + fieldseason_rain +
                            (1|p|SQUARE) +
                            ar(time = YRnm, gr = REP_ID)) + 
                       set_rescor(FALSE), data = mod_data, prior = mod_pr_red2,
                     save_pars = save_pars(all = TRUE, latent = TRUE), 
                     file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_WHW_PH_HAB",
                     cores = 4, iter = 4000, control = list(adapt_delta = 0.999))
summary(red_mod_2_whw)
plot(red_mod_2_whw, ask = FALSE)
pp_check(red_mod_2_whw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_whw, nsamples = 50, resp = "PH")

# plot(conditional_effects(red_mod_2_whw))
# check residuals against habitat
mod_resid <- resid(red_mod_2_whw)
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
ggsave("Residuals by habitat whole weighted reduced model 2.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")


loo_compare(full_mod_whw, red_mod_1_whw, red_mod_2_whw,
            criterion = "kfold")
# elpd_diff se_diff
# full_mod_whw    0.0       0.0  
# red_mod_1_whw -17.8      34.6  
# red_mod_2_whw -32.0      34.8

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

# Multivariate with N as a 2D spline
# Habitat interaction
# run model - full unweighted Ell R
full_mod_whuw <- update(full_mod_whw, newdata = mod_data,
                        control = list(adapt_delta = 0.95),
                        cores = 4, iter = 4000,
                        file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_WHUW_spl_PH_spl_HAB")

summary(full_mod_whuw)
plot(full_mod_whuw, ask = FALSE)
pp_check(full_mod_whuw, nsamples = 50, resp = "Ell")
pp_check(full_mod_whuw, nsamples = 50, resp = "PH")

# plot(conditional_effects(full_mod_whuw))
# check residuals against habitat
mod_resid <- resid(full_mod_whuw)
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
ggsave("Residuals by habitat whole unweighted full model.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")



# No pH ~ N spline
red_mod_1_whuw <- update(red_mod_1_whw, newdata = mod_data,
                         cores = 4, iter = 4000,
                         control = list(adapt_delta = 0.999),
                         file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_WHUW_spl_PH_HAB")

summary(red_mod_1_whuw)
plot(red_mod_1_whuw, ask = FALSE)
pp_check(red_mod_1_whuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_whuw, nsamples = 50, resp = "PH")

# plot(conditional_effects(red_mod_1_whuw))
# check residuals against habitat
mod_resid <- resid(red_mod_1_whuw)
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
ggsave("Residuals by habitat whole unweighted reduced model 1.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")

# No Ell ~ N spline
red_mod_2_whuw <- update(red_mod_2_whw, newdata = mod_data,
                         cores = 4, iter = 4000,
                         control = list(adapt_delta = 0.999),
                         file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_WHUW_PH_HAB")

summary(red_mod_2_whuw)
plot(red_mod_2_whuw, ask = FALSE)
pp_check(red_mod_2_whuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_whuw, nsamples = 50, resp = "PH")

# plot(conditional_effects(red_mod_2_whuw))
mod_resid <- resid(red_mod_2_whuw)
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
ggsave("Residuals by habitat whole unweighted reduced model 2.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")

loo_compare(full_mod_whuw, red_mod_1_whuw, red_mod_2_whuw,
            criterion = "kfold")
# elpd_diff se_diff
# red_mod_2_whuw   0.0       0.0  
# full_mod_whuw  -13.9      38.6  
# red_mod_1_whuw -17.0      17.6 

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
full_mod_smw <- update(full_mod_whw, newdata = mod_data,
                       cores = 4, iter = 5000,
                       control = list(adapt_delta = 0.95),
                       file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_SMW_spl_PH_spl_HAB")

summary(full_mod_smw)
plot(full_mod_smw, ask = FALSE)
pp_check(full_mod_smw, nsamples = 50, resp = "Ell")
pp_check(full_mod_smw, nsamples = 50, resp = "PH")

# plot(conditional_effects(full_mod_smw))
mod_resid <- resid(full_mod_smw)
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


# No pH ~ N spline
red_mod_1_smw <- update(red_mod_1_whw, newdata = mod_data,
                        cores = 4, iter = 4000,
                        control = list(adapt_delta = 0.99),
                        file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_SMW_spl_PH_HAB")

summary(red_mod_1_smw)
plot(red_mod_1_smw, ask = FALSE)
pp_check(red_mod_1_smw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_smw, nsamples = 50, resp = "PH")

# plot(conditional_effects(red_mod_1_smw))
mod_resid <- resid(red_mod_1_smw)
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
ggsave("Residuals by habitat small weighted reduced model 1.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")

# No Ell ~ N spline
red_mod_2_smw <- update(red_mod_2_whw, newdata = mod_data,
                        cores = 4, iter = 4000,
                        control = list(adapt_delta = 0.999),
                        file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_SMW_PH_HAB")

summary(red_mod_2_smw)
plot(red_mod_2_smw, ask = FALSE)
pp_check(red_mod_2_smw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_smw, nsamples = 50, resp = "PH")

# plot(conditional_effects(red_mod_2_smw))
mod_resid <- resid(red_mod_2_smw)
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
ggsave("Residuals by habitat small weighted reduced model 2.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")

loo_compare(full_mod_smw, red_mod_1_smw, red_mod_2_smw,
            criterion = "kfold")
# elpd_diff se_diff
# full_mod_smw     0.0       0.0 
# red_mod_2_smw  -93.3      33.4 
# red_mod_1_smw -113.1      32.6 

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
full_mod_smuw <- update(full_mod_whw, newdata = mod_data,
                        cores = 4, iter = 4000,
                        control = list(adapt_delta = 0.95),
                        file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_SMUW_spl_PH_spl_HAB")

summary(full_mod_smuw)
plot(full_mod_smuw, ask = FALSE)
pp_check(full_mod_smuw, nsamples = 50, resp = "Ell")
pp_check(full_mod_smuw, nsamples = 50, resp = "PH")

# plot(conditional_effects(full_mod_smuw))
mod_resid <- resid(full_mod_smuw)
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
ggsave("Residuals by habitat small unweighted full model.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")


# No pH ~ N spline
red_mod_1_smuw <- update(red_mod_1_whw, newdata = mod_data,
                         cores = 4, iter = 4000,
                         control = list(adapt_delta = 0.99),
                         file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_SMUW_spl_PH_HAB")

summary(red_mod_1_smuw)
plot(red_mod_1_smuw, ask = FALSE)
pp_check(red_mod_1_smuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_smuw, nsamples = 50, resp = "PH")

# plot(conditional_effects(red_mod_1_smuw))
mod_resid <- resid(red_mod_1_smuw)
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
ggsave("Residuals by habitat small unweighted reduced model 1.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")


# No Ell ~ N spline
red_mod_2_smuw <- update(red_mod_2_whw, newdata = mod_data,
                         cores = 4, iter = 4000,
                         control = list(adapt_delta = 0.999),
                         file = "Outputs/Models/Difference/Multivariate_measerrorXYU/ELL_SMUW_PH_HAB")

summary(red_mod_2_smuw)
plot(red_mod_2_smuw, ask = FALSE)
pp_check(red_mod_2_smuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_smuw, nsamples = 50, resp = "PH")

# plot(conditional_effects(red_mod_2_smuw))

mod_resid <- resid(red_mod_2_smuw)
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
ggsave("Residuals by habitat small unweighted reduced model 2.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU",
       width = 12, height = 18, units = "cm")

loo_compare(full_mod_smuw, red_mod_1_smuw, red_mod_2_smuw,
            criterion = "kfold")
# elpd_diff se_diff
# full_mod_smuw    0.0       0.0  
# red_mod_2_smuw -29.4      37.1  
# red_mod_1_smuw -55.5      38.5



# Summary plots ####
# Plot for comparing initial pH effects on pH change
nd <- 
  tibble(Year1_pH = seq(from = -2.25, to = 3.6, length.out = 30) %>% 
           rep(., times = 6),
         Ndep = rep(c(-1,0,1), each = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 90),
         Sdep = 0,
         fieldseason_rain = 0,
         PH_SE = 0.3,
         YRnm = 1,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
)) %>%
  mutate(Year1_pH = 1.2978*Year1_pH + 5.598,
         N = round(0.02433*N + 0.0743,2))

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
  mutate(N = as.character(N)) %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)) %>%
                filter(Management == "Low intensity"),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Initial pH", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("Initial pH and N effect on pH measerror under low intensity management unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 16, height = 12, units = "cm")

plot_dat %>%
  mutate(N = cut(scale(N), c(-10,-0.5,0.5,10), labels = c("0.05","0.07","0.1"))) %>%
  filter(!is.na(N)) %>%
  ggplot() +
  geom_point(aes(x = Year1_pH, y = PH, colour = N),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "Initial pH", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Initial pH and N effect on pH with data measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 16, height = 18, units = "cm")

# Plot for comparing Sdep and rainfall effects on pH change
nd <- 
  tibble(Sdep = seq(from = -5.5, to = 1.75, length.out = 30) %>% 
           rep(., times = 6),
         fieldseason_rain = rep(c(-1,0,1), each = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 90),
         Year1_pH = 0,
         N = 0,
         PH_SE = 0.3,
         YRnm = 1,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
)) %>%
  mutate(
    Sdep = 3.7523*Sdep - 4.8879,
    fieldseason_rain = round(30.436*fieldseason_rain + 7.914761)
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
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, fieldseason_rain = as.character(fieldseason_rain)),
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = fieldseason_rain, color = fieldseason_rain,
                  group = fieldseason_rain),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9, name = "Field season rain",
                         aesthetics = c("colour","fill"))
ggsave("Sdep and rain effect on pH measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 16, height = 18, units = "cm")

plot_dat %>%
  mutate(fieldseason_rain = cut(fieldseason_rain, c(-Inf,-15,23,Inf), labels = c("-23","8","38"))) %>%
  filter(!is.na(fieldseason_rain)) %>%
  ggplot() +
  geom_point(aes(x = Sdep, y = PH, colour = fieldseason_rain),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, fieldseason_rain = as.character(fieldseason_rain)),
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = fieldseason_rain, color = fieldseason_rain,
                  group = fieldseason_rain),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9, name = "Field season rain",
                         aesthetics = c("colour","fill"))
ggsave("Sdep and rain effect on pH with data measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 16, height = 18, units = "cm")

# Plot for comparing initial pH effects on Ellenberg R change
nd <- 
  tibble(Year1_pH = seq(from = -2.25, to = 3.6, length.out = 30) %>% 
           rep(., times = 6),
         N = rep(c(-1,0,1), each = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 90),
         Sdep = 0,
         PH = 0,
         ELL_SE = 0.4,
         fieldseason_rain = 0,
         YRnm = 1,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
)) %>%
  mutate(Year1_pH = 1.2978*Year1_pH + 5.598,
         N = round(0.02433*N + 0.0743,2))

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
  mutate(N = as.character(N)) %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)) %>%
                filter(Management == "Low intensity"),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_wrap(~Response) +
  labs(x = "Initial pH", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("Initial pH and N effect on Ellenberg R measerror low intensity management.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 16, height = 12, units = "cm")

plot_dat %>%
  mutate(N = cut(scale(N), c(-10,-0.5,0.5,10), labels = c("0.05","0.07","0.1"))) %>%
  filter(!is.na(N)) %>%
  ggplot() +
  geom_point(aes(x = Year1_pH, y = value, colour = N),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, N = as.character(N)),
              aes(x = Year1_pH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = N, color = N,
                  group = N),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "Initial pH", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Initial pH and N effect on Ellenberg R with data measerror unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 16, height = 18, units = "cm")

# Plot for comparing pH change on Ellenberg R change
nd <- 
  tibble(PH = seq(from = -3, to = 3.5, length.out = 40) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 40),
         Year1_pH = 0,
         N = 0,
         ELL_SE = 0.4,
         YRnm = 1,
         REP_ID = 1:80)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
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
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 20, height = 10, units = "cm")

plot_dat %>%
  filter(!is.na(N)) %>%
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
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 20, height = 10, units = "cm")

# ~~ no splines ####
# Plot for comparing pH change on Ellenberg R change
nd <- 
  tibble(PH = seq(from = -3, to = 3.5, length.out = 40) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 40),
         ELL_SE = 0.4,
         YRnm = 1,
         REP_ID = 1:80)

f <- do.call(rbind, list(
  fitted(red_mod_2_whw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(red_mod_2_whuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(red_mod_2_smw, newdata = nd, re_formula = NA, resp = "Ell") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(red_mod_2_smuw, newdata = nd, re_formula = NA, resp = "Ell") %>%
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
ggsave("pH change effect on Ellenberg R no spline multivariate model measerror.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 20, height = 10, units = "cm")

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_point(aes(x = PH, y = value, colour = Management),
             alpha = 1/3) +
  geom_smooth(data = f,
              aes(x = PH,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management),
              stat = "identity", 
              alpha = 1/2, size = 1/2) +
  facet_grid(~Response) +
  labs(x = "pH change", y = "Ellenberg R change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("pH change effect on Ellenberg R with data no spline model measerror.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 20, height = 10, units = "cm")

# Plot for comparing initial pH effects on pH change
nd <- 
  tibble(Sdep = seq(from = -5.5, to = 1.75, length.out = 30) %>% 
           rep(., times = 6),
         fieldseason_rain = rep(c(-1,0,1), each = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 90),
         YRnm = 1,
         PH_SE = 0.3,
         REP_ID = 1:180)

f <- do.call(rbind, list(
  fitted(red_mod_2_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(red_mod_2_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(red_mod_2_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(red_mod_2_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Sdep, fieldseason_rain) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"),
         Sdep =scale(Sdep),
         fieldseason_rain = scale(fieldseason_rain))

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, fieldseason_rain = as.character(fieldseason_rain)),
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = fieldseason_rain, color = fieldseason_rain,
                  group = fieldseason_rain),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Sdep and rain effect on pH no spline model measerror.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 16, height = 18, units = "cm")

plot_dat %>%
  mutate(fieldseason_rain = cut(scale(fieldseason_rain), c(-10,-0.5,0.5,10), labels = c("-1","0","1"))) %>%
  filter(!is.na(fieldseason_rain)) %>%
  ggplot() +
  geom_point(aes(x = Sdep, y = PH, colour = fieldseason_rain),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = mutate(f, fieldseason_rain = as.character(fieldseason_rain)),
              aes(x = Sdep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = fieldseason_rain, color = fieldseason_rain,
                  group = fieldseason_rain),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(Response~Management) +
  labs(x = "S deposition", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_viridis_d(option="plasma", end = 0.9) +
  scale_fill_viridis_d(option="plasma", end = 0.9)
ggsave("Sdep and rain effect on pH no spline model with data measerror.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 16, height = 18, units = "cm")

# plot for looking at fieldseason rain impact
nd <- 
  tibble(fieldseason_rain = seq(from = -3, to = 6.5, length.out = 30) %>% 
           rep(., times = 2),
         Management = rep(0:1, each = 30),
         Sdep = 0,
         YRnm = 1,
         PH_SE = 0.3,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(red_mod_2_whw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted full"),
  fitted(red_mod_2_whuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted full"),
  fitted(red_mod_2_smw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Weighted small"),
  fitted(red_mod_2_smuw, newdata = nd, re_formula = NA, resp = "PH") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Management = ifelse(Management == 1, "High intensity", "Low intensity"),
           Response = "Unweighted small")
))

plot_dat <- ELL_pH %>% 
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         PH,Management, REP_ID, Time_period,
         Sdep, fieldseason_rain) %>%
  filter(!is.na(Management)) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"),
         Management = ifelse(Management == "High", "High intensity", "Low intensity"),
         Sdep =scale(Sdep),
         fieldseason_rain = scale(fieldseason_rain))

plot_dat %>%
  ggplot() +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = fieldseason_rain,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity", 
              alpha = 1/4, size = 1/2) +
  facet_grid(~Response) +
  labs(x = "Field season rainfall difference", y = "pH change") +
  scale_x_continuous(expand = c(0,0)) 
ggsave("Rain effect on pH no spline model measerror.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 20, height = 12, units = "cm")

plot_dat %>%
  filter(!is.na(fieldseason_rain)) %>%
  ggplot() +
  geom_point(aes(x = fieldseason_rain, y = PH, colour = Management),
             alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = f,
              aes(x = fieldseason_rain,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = Management, color = Management,
                  group = Management),
              stat = "identity", 
              alpha = 1/3, size = 1/2) +
  facet_grid(~Response) +
  labs(x = "Field season rainfall difference", y = "pH change") +
  scale_x_continuous(expand = c(0,0))
ggsave("Rain effect on pH no spline model with data measerror.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 16, height = 18, units = "cm")

# Plot for comparing impact of Ndep upon N:C
nd <- 
  tibble(Ndep = seq(from = -1.8, to = 3.2, length.out = 60),
         YRnm = 1,
         REP_ID = 1:60)

f <- do.call(rbind, list(
  fitted(full_mod_whw, newdata = nd, re_formula = NA, resp = "N") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Response = "Weighted full"),
  fitted(full_mod_whuw, newdata = nd, re_formula = NA, resp = "N") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Response = "Unweighted full"),
  fitted(full_mod_smw, newdata = nd, re_formula = NA, resp = "N") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Response = "Weighted small"),
  fitted(full_mod_smuw, newdata = nd, re_formula = NA, resp = "N") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(Response = "Unweighted small")
))

f <- f %>%
  mutate(Estimate = Estimate*0.02436 + 0.0741,
         Q2.5 = Q2.5*0.02436 + 0.0741,
         Q97.5 = Q97.5*0.02436 + 0.0741,
         Ndep = Ndep*67.04 + 143.37)

plot_dat <- ELL_pH %>% 
  filter(!is.na(Management)) %>%
  select(WH_R_W,WH_R_UW,SM_R_W,SM_R_UW,
         N = NC_RATIO, Ndep, REP_ID, Time_period) %>%
  pivot_longer(contains("_R_"), names_to = "Response") %>%
  mutate(Response = recode(Response,
                           "WH_R_W" = "Weighted full",
                           "WH_R_UW" = "Unweighted full",
                           "SM_R_W" = "Weighted small",
                           "SM_R_UW" = "Unweighted small"))

plot_dat %>%
  filter(Response == "Unweighted small") %>%
  ggplot() +
  # geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = filter(f, Response == "Unweighted small"),
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity", 
              alpha = 1/4, size = 1/2,
              fill = "#CC79A7", colour = "#CC79A7") +
  labs(x = "N deposition", y = "N:C ratio") +
  scale_x_continuous(expand = c(0,0))
ggsave("Ndep effect on NC full model unweighted small unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 12, height = 12, units = "cm")

plot_dat %>%
  filter(Response == "Unweighted small") %>%
  ggplot() +
  geom_point(aes(x = Ndep, y = N),
             alpha = 0.5, colour = "#CC79A7") +
  # geom_hline(yintercept = 0, colour = "gray") +
  geom_smooth(data = filter(f, Response == "Unweighted small"),
              aes(x = Ndep,
                  y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity", 
              alpha = 1/3, size = 1/2,
              fill = "black", colour = "black") +
  labs(x = "N deposition", y = "N:C ratio") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits =c(0,0.16), expand = c(0,0))
ggsave("Ndep effect on NC full model unweighted small with data unscaled.png",
       path = "Outputs/Models/Difference/Multivariate_measerrorXYU", 
       width = 12, height = 12, units = "cm")

