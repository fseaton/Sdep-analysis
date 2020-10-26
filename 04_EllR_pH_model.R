# Ellenberg R and pH model
library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())

# taking PH_long and X_Ell from 01b script
ph_ell_long <- full_join(X_Ell, PH_long) %>%
  filter(!Year %in% c(1990,2016)) %>%
  mutate(SQNUM = sapply(strsplit(REP_ID, "X"), "[",1),
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
  left_join(select(CS_plot_atdep, REP_ID, Year, Ndep, Sdep) %>%
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
  left_join(select(CS_plot_atdep, REP_ID, Year, Ndep, Sdep) %>%
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
               ar(time = YRnm, gr = REP_ID),
             family = "student") +
            bf(PH  ~ Management*Sdep + fieldseason_rain + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), family = "student") + 
            bf(N ~ Ndep + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID)) +
            set_rescor(FALSE), data = mod_data)


mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.5), class = "b", resp = "N"),
            prior(student_t(3, 0, 2.5), class = "sds", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(student_t(3,0,1), class = "Intercept", resp = "N"),
            prior(gamma(4,1), class = "nu", resp = "Ell"),
            prior(gamma(4,1), class = "nu", resp = "PH"),
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
                      ar(time = YRnm, gr = REP_ID),
                    family = "student") +
                   bf(PH  ~ Management*Sdep + fieldseason_rain  + s(Year1_pH, N, Management) +
                        (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID), family = "student") + 
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
full_mod_whw <- brm(bf(Ell  ~ Management*PH + s(Year1_pH, N, Management) + 
                         (1|SQUARE) +
                         ar(time = YRnm, gr = REP_ID),
                       family = "student") +
                      bf(PH  ~ Management*Sdep + fieldseason_rain  + s(Year1_pH, N, Management) +
                           (1|SQUARE) +
                           ar(time = YRnm, gr = REP_ID), family = "student") + 
                      bf(N ~ Ndep + (1|SQUARE) +
                           ar(time = YRnm, gr = REP_ID)) +
                      set_rescor(FALSE), data = mod_data, prior = mod_pr,
                    save_pars = save_pars(all = TRUE, latent = TRUE), 
                    file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHW_spl_PH_spl_N_HAB",
                    cores = 4, iter = 5000, control = list(adapt_delta = 0.95))
summary(full_mod_whw)
plot(full_mod_whw, ask = FALSE)
pp_check(full_mod_whw, nsamples = 50, resp = "Ell")
pp_check(full_mod_whw, nsamples = 50, resp = "PH")
pp_check(full_mod_whw, nsamples = 50, resp = "N")

# plot(conditional_effects(full_mod_whw))

full_mod_whw <- add_criterion(full_mod_whw, "loo",  moment_match = TRUE, reloo = TRUE,
                              file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHW_spl_PH_spl_N_HAB")

# No pH ~ N spline
mod_pr_red1 <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
                 prior(normal(0,0.5), class = "b", resp = "PH"),
                 prior(normal(0,0.5), class = "b", resp = "N"),
                 prior(student_t(3, 0, 2.5), class = "sds", resp = "Ell"),
                 prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
                 prior(normal(0,0.25), class = "Intercept", resp = "PH"),
                 prior(student_t(3,0,1), class = "Intercept", resp = "N"),
                 prior(gamma(4,1), class = "nu", resp = "Ell"),
                 prior(gamma(4,1), class = "nu", resp = "PH"),
                 prior(normal(0,0.2), class = "ar", resp = "Ell"),
                 prior(normal(0,0.2), class = "ar", resp = "PH"),
                 prior(normal(0,0.2), class = "ar", resp = "N"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sd", resp = "N"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
                 prior(student_t(3,0,0.5), class = "sigma", resp = "N"))

red_mod_1_whw <- brm(bf(Ell  ~ Management*PH + s(Year1_pH, N, Management) + 
                          (1|SQUARE) +
                          ar(time = YRnm, gr = REP_ID),
                        family = "student") +
                       bf(PH  ~ Management*Sdep + fieldseason_rain +
                            (1|SQUARE) +
                            ar(time = YRnm, gr = REP_ID), family = "student") + 
                       bf(N ~ Ndep + (1|SQUARE) +
                            ar(time = YRnm, gr = REP_ID)) +
                       set_rescor(FALSE), data = mod_data, prior = mod_pr_red1,
                     save_pars = save_pars(all = TRUE, latent = TRUE), 
                     file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHW_spl_PH_N_HAB",
                     cores = 4, iter = 5000, control = list(adapt_delta = 0.95))
summary(red_mod_1_whw)
plot(red_mod_1_whw, ask = FALSE)
pp_check(red_mod_1_whw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_whw, nsamples = 50, resp = "PH")
pp_check(red_mod_1_whw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_1_whw))

red_mod_1_whw <- add_criterion(red_mod_1_whw, "loo",  moment_match = TRUE, reloo = TRUE,
                               file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHW_spl_PH_N_HAB")

# No Ell ~ N spline
mod_pr_red2 <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
                 prior(normal(0,0.5), class = "b", resp = "PH"),
                 prior(normal(0,0.5), class = "b", resp = "N"),
                 prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
                 prior(normal(0,0.25), class = "Intercept", resp = "PH"),
                 prior(student_t(3,0,1), class = "Intercept", resp = "N"),
                 prior(gamma(4,1), class = "nu", resp = "Ell"),
                 prior(gamma(4,1), class = "nu", resp = "PH"),
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
                          ar(time = YRnm, gr = REP_ID),
                        family = "student") +
                       bf(PH  ~ Management*Sdep + fieldseason_rain +
                            (1|SQUARE) +
                            ar(time = YRnm, gr = REP_ID), family = "student") + 
                       bf(N ~ Ndep + (1|SQUARE) +
                            ar(time = YRnm, gr = REP_ID)) +
                       set_rescor(FALSE), data = mod_data, prior = mod_pr_red2,
                     save_pars = save_pars(all = TRUE, latent = TRUE), 
                     file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHW_PH_N_HAB",
                     cores = 4, iter = 5000, control = list(adapt_delta = 0.95))
summary(red_mod_2_whw)
plot(red_mod_2_whw, ask = FALSE)
pp_check(red_mod_2_whw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_whw, nsamples = 50, resp = "PH")
pp_check(red_mod_2_whw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_2_whw))

red_mod_2_whw <- add_criterion(red_mod_2_whw, "loo",  moment_match = TRUE, reloo = TRUE, 
                               file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHW_PH_N_HAB")

loo_compare(full_mod_whw, red_mod_1_whw, red_mod_2_whw)

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
# run model - full weighted Ell R
full_mod_whuw <- update(full_mod_whw, newdata = mod_data,
                        control = list(adapt_delta = 0.95),
                        cores = 4, iter = 5000,
                        file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHUW_spl_PH_spl_N_HAB")

summary(full_mod_whuw)
plot(full_mod_whuw, ask = FALSE)
pp_check(full_mod_whuw, nsamples = 50, resp = "Ell")
pp_check(full_mod_whuw, nsamples = 50, resp = "PH")
pp_check(full_mod_whuw, nsamples = 50, resp = "N")

# plot(conditional_effects(full_mod_whuw))

full_mod_whuw <- add_criterion(full_mod_whuw, "loo",  moment_match = TRUE, reloo = TRUE,
                               file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHUW_spl_PH_spl_N_HAB")

# No pH ~ N spline
red_mod_1_whuw <- update(red_mod_1_whw, newdata = mod_data,
                         cores = 4, iter = 5000,
                         control = list(adapt_delta = 0.95),
                         file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHUW_spl_PH_N_HAB")

summary(red_mod_1_whuw)
plot(red_mod_1_whuw, ask = FALSE)
pp_check(red_mod_1_whuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_whuw, nsamples = 50, resp = "PH")
pp_check(red_mod_1_whuw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_1_whuw))

red_mod_1_whuw <- add_criterion(red_mod_1_whuw, "loo",  moment_match = TRUE, reloo = TRUE,
                                file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHUW_spl_PH_N_HAB")
# No Ell ~ N spline
red_mod_2_whuw <- update(red_mod_2_whw, newdata = mod_data,
                         cores = 4, iter = 5000,
                         control = list(adapt_delta = 0.95),
                         file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHUW_PH_N_HAB")

summary(red_mod_2_whuw)
plot(red_mod_2_whuw, ask = FALSE)
pp_check(red_mod_2_whuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_whuw, nsamples = 50, resp = "PH")
pp_check(red_mod_2_whuw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_2_whuw))

red_mod_2_whuw <- add_criterion(red_mod_2_whuw, "loo",  moment_match = TRUE, reloo = TRUE,
                                file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_WHUW_PH_N_HAB")

loo_compare(full_mod_whuw, red_mod_1_whuw, red_mod_2_whuw)

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
# run model - full weighted Ell R
full_mod_smw <- update(full_mod_whw, newdata = mod_data,
                       cores = 4, iter = 5000,
                       control = list(adapt_delta = 0.95),
                       file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMW_spl_PH_spl_N_HAB")

summary(full_mod_smw)
plot(full_mod_smw, ask = FALSE)
pp_check(full_mod_smw, nsamples = 50, resp = "Ell")
pp_check(full_mod_smw, nsamples = 50, resp = "PH")
pp_check(full_mod_smw, nsamples = 50, resp = "N")

# plot(conditional_effects(full_mod_smw))

full_mod_smw <- add_criterion(full_mod_smw, "loo",  moment_match = TRUE, reloo = TRUE,
                              file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMW_spl_PH_spl_N_HAB")

# No pH ~ N spline
red_mod_1_smw <- update(red_mod_1_whw, newdata = mod_data,
                        cores = 4, iter = 5000,
                        control = list(adapt_delta = 0.95),
                        file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMW_spl_PH_N_HAB")

summary(red_mod_1_smw)
plot(red_mod_1_smw, ask = FALSE)
pp_check(red_mod_1_smw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_smw, nsamples = 50, resp = "PH")
pp_check(red_mod_1_smw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_1_smw))

red_mod_1_smw <- add_criterion(red_mod_1_smw, "loo",  moment_match = TRUE, reloo = TRUE,
                               file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMW_spl_PH_N_HAB")
# No Ell ~ N spline
red_mod_2_smw <- update(red_mod_2_whw, newdata = mod_data,
                        cores = 4, iter = 5000,
                        control = list(adapt_delta = 0.95),
                        file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMW_PH_N_HAB")

summary(red_mod_2_smw)
plot(red_mod_2_smw, ask = FALSE)
pp_check(red_mod_2_smw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_smw, nsamples = 50, resp = "PH")
pp_check(red_mod_2_smw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_2_smw))

red_mod_2_smw <- add_criterion(red_mod_2_smw, "loo",  moment_match = TRUE, reloo = TRUE,
                               file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMW_PH_N_HAB")

loo_compare(full_mod_smw, red_mod_1_smw, red_mod_2_smw)

# ~~~ 4m2 unweighted ####
mod_data <- ELL_pH %>%
  filter(!is.na(Management)) %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW) %>%
  select(REP_ID, SQUARE, YRnm, Management, Ell, Sdep, Ndep,
         fieldseason_rain, Year2_pH, PH, 
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
# run model - full weighted Ell R
full_mod_smuw <- update(full_mod_whw, newdata = mod_data,
                        cores = 4, iter = 5000,
                        control = list(adapt_delta = 0.95),
                        file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMUW_spl_PH_spl_N_HAB")

summary(full_mod_smuw)
plot(full_mod_smuw, ask = FALSE)
pp_check(full_mod_smuw, nsamples = 50, resp = "Ell")
pp_check(full_mod_smuw, nsamples = 50, resp = "PH")
pp_check(full_mod_smuw, nsamples = 50, resp = "N")

# plot(conditional_effects(full_mod_smuw))

full_mod_smuw <- add_criterion(full_mod_smuw, "loo",  moment_match = TRUE, reloo = TRUE,
                               file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMUW_spl_PH_spl_N_HAB")

# No pH ~ N spline
red_mod_1_smuw <- update(red_mod_1_whw, newdata = mod_data,
                         cores = 4, iter = 5000,
                         control = list(adapt_delta = 0.95),
                         file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMUW_spl_PH_N_HAB")

summary(red_mod_1_smuw)
plot(red_mod_1_smuw, ask = FALSE)
pp_check(red_mod_1_smuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_1_smuw, nsamples = 50, resp = "PH")
pp_check(red_mod_1_smuw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_1_smuw))

red_mod_1_smuw <- add_criterion(red_mod_1_smuw, "loo",  moment_match = TRUE, reloo = TRUE,
                                file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMUW_spl_PH_N_HAB")
# No Ell ~ N spline
red_mod_2_smuw <- update(red_mod_2_whw, newdata = mod_data,
                         cores = 4, iter = 5000,
                         control = list(adapt_delta = 0.95),
                         file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMUW_PH_N_HAB")

summary(red_mod_2_smuw)
plot(red_mod_2_smuw, ask = FALSE)
pp_check(red_mod_2_smuw, nsamples = 50, resp = "Ell")
pp_check(red_mod_2_smuw, nsamples = 50, resp = "PH")
pp_check(red_mod_2_smuw, nsamples = 50, resp = "N")

# plot(conditional_effects(red_mod_2_smuw))

red_mod_2_smuw <- add_criterion(red_mod_2_smuw, "loo",  moment_match = TRUE, reloo = TRUE,
                                file = "Outputs/Models/Difference/Multivariate_nomeaserror/ELL_SMUW_PH_N_HAB")

loo_compare(full_mod_smuw, red_mod_1_smuw, red_mod_2_smuw)


# ** Measurement error ####
# Multivariate with N as a 2D spline
# Cannot have mi notation within spline
get_prior(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Year1_pH, N) + (1|p|SQUARE) +
               ar(time = YRnm, gr = REP_ID),
             family = "student") +
            bf(PH | mi(PH_DIW_SE)  ~ Sdep + fieldseason_rain + (1|p|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), family = "student") + 
            bf(N ~ Ndep + C + (1|p|SQUARE) +
                 ar(time = YRnm, gr = REP_ID)) +
            set_rescor(FALSE), data = mod_data)


mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.5), class = "b", resp = "N"),
            prior(student_t(3, 0, 2.5), class = "sds", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(student_t(3,0,1), class = "Intercept", resp = "N"),
            prior(gamma(4,1), class = "nu", resp = "Ell"),
            prior(gamma(4,1), class = "nu", resp = "PH"),
            prior(normal(0,0.5), class = "ar", resp = "Ell"),
            prior(normal(0,0.5), class = "ar", resp = "PH"),
            prior(normal(0,0.5), class = "ar", resp = "N"),
            prior(student_t(3,0,0.5), class = "sd", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sd", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sd", resp = "N"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "Ell"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "PH"),
            prior(student_t(3,0,0.5), class = "sigma", resp = "N"),
            prior(lkj(2), class = "cor", group = "SQUARE"))
make_stancode(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Year1_pH, N) + (1|p|SQUARE) +
                   ar(time = YRnm, gr = REP_ID),
                 family = "student") +
                bf(PH | mi(PH_DIW_SE)  ~ Sdep + fieldseason_rain + (1|p|SQUARE) +
                     ar(time = YRnm, gr = REP_ID), family = "student") + 
                bf(N ~ Ndep + C + (1|p|SQUARE) +
                     ar(time = YRnm, gr = REP_ID)) +
                set_rescor(FALSE), data = mod_data, prior = mod_pr)

# prior simulation
prior_mod <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(Year1_pH, N) + (1|p|SQUARE) +
                      ar(time = YRnm, gr = REP_ID),
                    family = "student") +
                   bf(PH | mi(PH_DIW_SE)  ~ Sdep + fieldseason_rain + (1|p|SQUARE) +
                        ar(time = YRnm, gr = REP_ID), family = "student") + 
                   bf(N ~ Ndep + C + (1|p|SQUARE) +
                        ar(time = YRnm, gr = REP_ID)) +
                   set_rescor(FALSE), data = mod_data, prior = mod_pr,
                 sample_prior = "only", save_pars = save_pars(all = TRUE, latent = TRUE), 
                 cores = 4, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50, resp = "Ell")
pp_check(prior_mod, nsamples = 50, resp = "PH")
plot(conditional_effects(prior_mod))

# randomly pick 300 rows and run model to see what happens
mod_data2 <- mod_data[sample.int(nrow(mod_data),300),]
full_mod <- brm(bf(Ell | mi(ELL_SE) ~ mi(PH) + s(N, Year1_pH) + (1|p|SQUARE) +
                     ar(time = YRnm, gr = REP_ID),
                   family = "student") +
                  bf(PH | mi(PH_DIW_SE)  ~ Sdep + fieldseason_rain + (1|p|SQUARE) +
                       ar(time = YRnm, gr = REP_ID), family = "student") + 
                  bf(N ~ Ndep + C + (1|p|SQUARE) +
                       ar(time = YRnm, gr = REP_ID)) +
                  set_rescor(FALSE), data = mod_data2, prior = mod_pr,
                save_pars = save_pars(all = TRUE, latent = TRUE), 
                cores = 4, iter = 5000)
summary(full_mod)
pp_check(full_mod, nsamples = 50, resp = "Ell")
pp_check(full_mod, nsamples = 50, resp = "PH")
pp_check(full_mod, nsamples = 50, resp = "N")

plot(conditional_effects(full_mod))

# check where 6 divergent transitions are
mypars <- colnames(as.matrix(full_mod))
pairs(full_mod$fit, pars = grep("^b_|^sigma_|^nu_", mypars, value = TRUE))
# lower triangle so upping adapt delta should fix it


# Habitat interaction
mod_data <- mutate(mod_data, Management = ifelse(Management == "High",1,0))
get_prior(bf(Ell | mi(ELL_SE) ~ Management*mi(PH) + s(Year2_pH, N, Management) + 
               (1|SQUARE) +
               ar(time = YRnm, gr = REP_ID),
             family = "student") +
            bf(PH | mi(PH_DIW_SE)  ~ Management*Sdep + fieldseason_rain + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), family = "student") + 
            bf(N ~ Ndep + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID)) +
            set_rescor(FALSE), data = mod_data)


mod_pr <- c(prior(normal(0,0.5), class = "b", resp = "Ell"),
            prior(normal(0,0.5), class = "b", resp = "PH"),
            prior(normal(0,0.5), class = "b", resp = "N"),
            prior(student_t(3, 0, 2.5), class = "sds", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "Ell"),
            prior(normal(0,0.25), class = "Intercept", resp = "PH"),
            prior(student_t(3,0,1), class = "Intercept", resp = "N"),
            prior(gamma(4,1), class = "nu", resp = "Ell"),
            prior(gamma(4,1), class = "nu", resp = "PH"),
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
prior_mod <- brm(bf(Ell | mi(ELL_SE) ~ Management*mi(PH) + s(Year2_pH, N, Management) + 
                      (1|SQUARE) +
                      ar(time = YRnm, gr = REP_ID),
                    family = "student") +
                   bf(PH | mi(PH_DIW_SE)  ~ Management*Sdep + fieldseason_rain + (1|SQUARE) +
                        ar(time = YRnm, gr = REP_ID), family = "student") + 
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

# run model to see what happens
full_mod <- brm(bf(Ell | mi(ELL_SE) ~ Management*mi(PH) + s(Year2_pH, N, Management) + 
                     (1|SQUARE) +
                     ar(time = YRnm, gr = REP_ID),
                   family = "student") +
                  bf(PH | mi(PH_DIW_SE)  ~ Management*Sdep + fieldseason_rain + (1|SQUARE) +
                       ar(time = YRnm, gr = REP_ID), family = "student") + 
                  bf(N ~ Ndep + (1|SQUARE) +
                       ar(time = YRnm, gr = REP_ID)) +
                  set_rescor(FALSE), data = mod_data, prior = mod_pr,
                save_pars = save_pars(all = TRUE, latent = TRUE), 
                cores = 4, iter = 10000, control = list(adapt_delta = 0.99))
summary(full_mod)
pp_check(full_mod, nsamples = 50, resp = "Ell")
pp_check(full_mod, nsamples = 50, resp = "PH")
pp_check(full_mod, nsamples = 50, resp = "N")

plot(conditional_effects(full_mod))

ell_eutr_conditions <- data.frame(
  cond__ = c("Management = 0 & Year2_pH = 4.44",
             "Management = 0 & Year2_pH = 5.55",
             "Management = 0 & Year2_pH = 6.66",
             "Management = 1 & Year2_pH = 4.44",
             "Management = 1 & Year2_pH = 5.55",
             "Management = 1 & Year2_pH = 6.66")
)
plot(conditional_effects(full_mod, effects = "N",
                         conditions = ell_eutr_conditions))

# simulate data according to model and check if works 
set.seed(151020)
sim_data <- data.frame(SQUARE = gl(100,5)) %>%
  mutate(REP_ID = paste0(SQUARE,"X",1:500),    
         Improved = rbinom(500,1,0.007*as.numeric(SQUARE))) %>%
  mutate(pH = rnorm(500,mean = 4.5 + Improved + rep(rnorm(100,0,0.5),each = 5), 1)) %>%
  mutate(Ell = rnorm(500, mean = pH + rep(rnorm(100,0,1),each = 5), 1)) %>%
  mutate(logC_year1 = rnorm(500, rep(rnorm(100,0,0.2),each = 5), 1)) %>%
  mutate(logN_year1 = rnorm(500, logC_year1, 0.02)) %>%
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
                                      (pH>5.5)*-0.5*logN_year1,EllSE_year1)) %>%
  mutate(pH_year2 = pH + pH_diffYear2,
         pHSE_year2 = 0.15,
         Ell_year2 = Ell + Ell_diffYear2,
         EllSE_year2 = 0.2 + Improved*0.05,
         logN_year2 = logN_year1 + rnorm(500,0.5*(2+Ndep),0.2),
         logC_year2 = logC_year1 + rnorm(500,0,0.2),
         rain_year3 = rep(rnorm(100,0,1),each = 5)) %>%
  mutate(pH_diffYear3 = rstudent_t(500, 4,
                                   0.5*Improved+(1-Improved)*0.5*Sdep + 
                                     0.3*rain_year3, pHSE_year2)) %>%
  mutate(Ell_diffYear3 = rstudent_t(500, 4,
                                    pH_diffYear3*(1-Improved) + 
                                      Improved*0.5*pH_diffYear3 + 
                                      (pH_year2>5.5)*-0.5*logN_year2,EllSE_year2)) %>%
  select(SQUARE, REP_ID, Sdep, Ndep, Improved,
         rain_year12 = rain_year2,rain_year23 = rain_year3,
         C_year12 = logC_year1, C_year23 = logC_year2,
         N_year12 = logN_year1, N_year23 = logN_year2,
         PH_DIW_SE_year12 = pHSE_year1, PH_DIW_SE_year23 = pHSE_year2,
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


sim_mod <- brm(bf(Ell | mi(ELL_SE) ~ Improved*mi(PH) + s(Year1_pH, N, Improved) + 
                    (1|SQUARE) +
                    ar(time = YRnm, gr = REP_ID),
                  family = "student") +
                 bf(PH | mi(PH_DIW_SE)  ~ Improved*Sdep + rain + (1|SQUARE) +
                      ar(time = YRnm, gr = REP_ID), family = "student") + 
                 bf(N ~ Ndep + C + (1|SQUARE) +
                      ar(time = YRnm, gr = REP_ID)) +
                 set_rescor(FALSE), data = sim_data, prior = mod_pr,
               save_pars = save_pars(all = TRUE, latent = TRUE), 
               cores = 4, iter = 4000)
summary(sim_mod)
pp_check(sim_mod, nsamples = 50, resp = "Ell")
pp_check(sim_mod, nsamples = 50, resp = "PH")
pp_check(sim_mod, nsamples = 50, resp = "N")

plot(conditional_effects(sim_mod))



get_prior(bf(Ell | mi(ELL_SE) ~ Inter + b1*mi(PH) + b2*step(Year1_pH - phthreshold)*N,
             b1 + b2 + phthreshold ~ 1,Inter ~ (1|p|SQUARE),# + ar(time = YRnm, gr = REP_ID),
             nl = TRUE, family = "student") +
            bf(PH | mi(PH_DIW_SE)  ~ Sdep + mi(Moisture) + (1|p|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), family = "student") + 
            bf(Moisture | mi() ~ fieldseason_rain*LOI + (1|SQUARE) +
                 ar(time = Yrnm, gr = REP_ID)) +
            bf(N ~ Ndep + (1|p|SQUARE) +
                 ar(time = YRnm, gr = REP_ID)) +
            set_rescor(FALSE), data = mod_data)