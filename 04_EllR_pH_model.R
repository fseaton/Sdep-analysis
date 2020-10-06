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
                   data = mod_data, prior = pr, cores = 6, chains = 4,
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
              data = mod_data, prior = pr, cores = 6, chains = 4)
summary(mod_nl)
plot(mod_nl)
plot(conditional_effects(mod_pr_only))

pp_check(mod_pr_only)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# Ellenberg and pH difference models ####

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
# 39 missing 1998 but having 2007/1978, 23 missing 2007 but 
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
                       "0719" = 2019)) %>%
  left_join(PH_QA_diff) %>%
  left_join(BH_comb_nodupes) %>%
  left_join(select(CS_plot_atdep, REP_ID, Year, Ndep, Sdep)) %>%
  left_join(cs_rainfall_diff)


janitor::get_dupes(ELL_pH, Time_period, REP_ID) # no duplicates
mice::md.pattern(ELL_pH)
View(filter(ELL_pH, is.na(WH_R_UW) & !is.na(PH)))



# simple difference model
# model each Ellenberg R change as a function of pH change 
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W) %>%
  filter(!BH %in% c(4,3,2,5,20,21,14,22,13,18,19,17)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH) %>%
  na.omit()

# normal model response - to show not great
get_prior(Ell ~ PH + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
          data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(normal(0.4,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))
normal_mod <- brm(Ell ~ PH + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                  data = mod_data, prior = mod_pr,
                  cores = 6, iter = 5000)
pp_check(test_mod, nsamples = 50) + scale_x_continuous(limits =c(-5,5))

# prior checks
get_prior(Ell ~ PH + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
          family = "student", data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(gamma(2,0.1), class = "nu"),
            prior(normal(0.4,0.2), class = "ar"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))

# prior predictive check
prior_mod <- brm(Ell ~ PH + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                 family = "student", sample_prior = "only",
                 data = mod_data, prior = mod_pr, cores = 6, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50)
# underestimates peak at 0 on average but gets close

# ellenberg R model - weighted Ellenberg R for whole plot
ell_ph_diff <- brm(Ell ~ PH + (1|SQUARE) + ar(time = YRnm, gr = REP_ID),
                   data = mod_data, prior = mod_pr, family = "student",
                   cores = 6, iter = 5000, 
                   file = "Outputs/Models/Difference/EllR_WHW_PH")
summary(ell_ph_diff)
plot(ell_ph_diff, plot = FALSE)
pp_check(ell_ph_diff)

# unweighted ellenberg r whole plot
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW) %>%
  filter(!BH %in% c(4,3,2,5,20,21,14,22,13,18,19,17)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH) %>%
  na.omit() 
ell_ph_diffuw <- update(ell_ph_diff, newdata = mod_data, cores=6, iter = 5000,
                        file = "Outputs/Models/Difference/EllR_WHUW_PH",
                        control = list(adapt_delta = 0.95))
summary(ell_ph_diffuw)
plot(ell_ph_diffuw, plot = FALSE)
pp_check(ell_ph_diffuw)


# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W) %>%
  filter(!BH %in% c(4,3,2,5,20,21,14,22,13,18,19,17)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH) %>%
  na.omit() 
ell_ph_diffsmw <- update(ell_ph_diff, newdata = mod_data, cores=6, iter = 5000,
                         file = "Outputs/Models/Difference/EllR_SMW_PH",
                         control = list(adapt_delta = 0.95))
summary(ell_ph_diffsmw)
plot(ell_ph_diffsmw, plot = FALSE)
pp_check(ell_ph_diffsmw)

# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW) %>%
  filter(!BH %in% c(4,3,2,5,20,21,14,22,13,18,19,17)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH) %>%
  na.omit() 
ell_ph_diffsmuw <- update(ell_ph_diff, newdata = mod_data, cores=6, iter = 5000,
                          control = list(adapt_delta = 0.95),
                          file = "Outputs/Models/Difference/EllR_SMUW_PH")
summary(ell_ph_diffsmuw)
plot(ell_ph_diffsmuw, plot = FALSE)
pp_check(ell_ph_diffsmuw)


# Plot for comparing pH effects on different Ellenberg R scores
ph_fixef <- rbind(fixef(ell_ph_diff),
                  fixef(ell_ph_diffuw),
                  fixef(ell_ph_diffsmw),
                  fixef(ell_ph_diffsmuw)) %>%
  as.data.frame() %>%
  cbind(score = rep(c("Weighted full",
                      "Unweighted full",
                      "Weighted small",
                      "Unweighted small"), each = 2)) %>%
  cbind(Variable = rep(c("Intercept","pH"), 4)) %>%
  as.data.frame()

ggplot(filter(ph_fixef, Variable == "pH"), 
       aes(x = score, y = Estimate)) +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5)) +
  geom_linerange(aes(ymin = Estimate - Est.Error, 
                     ymax = Estimate + Est.Error),
                 size = 1.5) + 
  scale_y_continuous(limits =c(0,0.11), expand = c(0,0)) +
  labs(x = "", y = "Effect of change in pH upon Ellenberg R")
ggsave("pH effect on Ellenberg R score versions.png",
       path = "Outputs/Models/Difference", 
       width = 12, height = 12, units = "cm")

# ** pH diw vs cacl2 comparison ####
# only one time transition so no temporal autocorrelation included
mod_data <- ELL_pH %>%
  mutate(SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W) %>%
  filter(!BH %in% c(4,3,2,5,20,21,14,22,13,18,19,17)) %>%
  select(REP_ID, SQUARE, Ell, PH, PHC) %>%
  na.omit()

get_prior(Ell ~ PH + (1|SQUARE),
          data = mod_data)

mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(0,0.25), class = "Intercept"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"))

ell_ph_diff <- brm(Ell ~ PH + (1|SQUARE), family = "student",
                   data = mod_data, prior = mod_pr, cores = 6, iter = 5000, 
                   file = "Outputs/Models/Difference/EllR_WHW_PH_0719only",
                   save_all_pars = TRUE, control = list(adapt_delta = 0.95))
summary(ell_ph_diff)
plot(ell_ph_diff)
pp_check(ell_ph_diff)
ell_ph_diff <- add_criterion(ell_ph_diff, "loo", moment_match = TRUE, reloo = TRUE)

ell_phc_diff <- brm(Ell ~ PHC + (1|SQUARE), family = "student",
                    data = mod_data, prior = mod_pr, cores = 6, iter = 5000, 
                    file = "Outputs/Models/Difference/EllR_WHW_PHC_0719only",
                    save_all_pars = TRUE,  control = list(adapt_delta = 0.95))
summary(ell_phc_diff)
plot(ell_phc_diff)
pp_check(ell_phc_diff)
ell_phc_diff <- add_criterion(ell_phc_diff, "loo", moment_match = TRUE, reloo = TRUE,
                              file = "Outputs/Models/Difference/EllR_WHW_PHC_0719only")

loo_compare(ell_ph_diff, ell_phc_diff)
#               elpd_diff se_diff
# ell_phc_diff  0.0       0.0   
# ell_ph_diff  -0.5       1.2

# unweighted ellenberg r whole plot
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_UW) %>%
  filter(!BH %in% c(4,3,2,5,20,21,14,22,13,18,19,17)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, PHC) %>%
  na.omit()
ell_ph_diffuw <- update(ell_ph_diff, newdata = mod_data, cores=6, iter = 5000,
                        file = "Outputs/Models/Difference/EllR_WHUW_PH_0719only")
summary(ell_ph_diffuw)
plot(ell_ph_diffuw)
pp_check(ell_ph_diffuw)
ell_ph_diffuw <- add_criterion(ell_ph_diffuw, "loo", moment_match = TRUE, reloo = TRUE,
                               file = "Outputs/Models/Difference/EllR_WHUW_PH_0719only")

ell_phc_diffuw <- update(ell_phc_diff, newdata = mod_data, cores=6, iter = 5000,
                         file = "Outputs/Models/Difference/EllR_WHUW_PHC_0719only")
summary(ell_phc_diffuw)
plot(ell_phc_diffuw)
pp_check(ell_phc_diffuw)
ell_phc_diffuw <- add_criterion(ell_phc_diffuw, "loo", moment_match = TRUE, reloo = TRUE,
                                file = "Outputs/Models/Difference/EllR_WHUW_PHC_0719only")

loo_compare(ell_ph_diffuw, ell_phc_diffuw)
# elpd_diff se_diff
# ell_ph_diffuw   0.0       0.0   
# ell_phc_diffuw -0.4       0.5  

# weighted ellenberg r small plot
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_W) %>%
  filter(!BH %in% c(4,3,2,5,20,21,14,22,13,18,19,17)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, PHC) %>%
  na.omit()
ell_ph_diffsmw <- update(ell_ph_diff, newdata = mod_data, cores=6, iter = 5000,
                         file = "Outputs/Models/Difference/EllR_SMW_PH_0719only")
summary(ell_ph_diffsmw)
plot(ell_ph_diffsmw)
pp_check(ell_ph_diffsmw)
ell_ph_diffsmw <- add_criterion(ell_ph_diffsmw, "loo", moment_match = TRUE, reloo = TRUE,
                                file = "Outputs/Models/Difference/EllR_SMW_PH_0719only")

ell_phc_diffsmw <- update(ell_phc_diff, newdata = mod_data, cores=6, iter = 5000,
                          file = "Outputs/Models/Difference/EllR_SMW_PHC_0719only")
summary(ell_phc_diffsmw)
plot(ell_phc_diffsmw)
pp_check(ell_phc_diffsmw)
ell_phc_diffsmw <- add_criterion(ell_phc_diffsmw, "loo", moment_match = TRUE, reloo = TRUE,
                                 file = "Outputs/Models/Difference/EllR_SMW_PHC_0719only")

loo_compare(ell_ph_diffsmw, ell_phc_diffsmw)
#                 elpd_diff se_diff
# ell_ph_diffsmw   0.0       0.0   
# ell_phc_diffsmw -0.8       1.3

# unweighted ellenberg r small plot
mod_data <- ELL_pH %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = SM_R_UW) %>%
  filter(!BH %in% c(4,3,2,5,20,21,14,22,13,18,19,17)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, PH, PHC) %>%
  na.omit()
ell_ph_diffsmuw <- update(ell_ph_diff, newdata = mod_data, cores=6, iter = 5000,
                          file = "Outputs/Models/Difference/EllR_SMUW_PH_0719only")
summary(ell_ph_diffsmuw)
plot(ell_ph_diffsmuw)
pp_check(ell_ph_diffsmuw)
ell_ph_diffsmuw <- add_criterion(ell_ph_diffsmuw, "loo", moment_match = TRUE, reloo = TRUE,
                                 file = "Outputs/Models/Difference/EllR_SMUW_PH_0719only")

ell_phc_diffsmuw <- update(ell_phc_diff, newdata = mod_data, cores=6, iter = 5000,
                           file = "Outputs/Models/Difference/EllR_SMUW_PHC_0719only")
summary(ell_phc_diffsmuw)
plot(ell_phc_diffsmuw)
pp_check(ell_phc_diffsmuw)
ell_phc_diffsmuw <- add_criterion(ell_phc_diffsmuw, "loo", moment_match = TRUE, reloo = TRUE,
                                  file = "Outputs/Models/Difference/EllR_SMUW_PHC_0719only")

loo_compare(ell_ph_diffsmuw, ell_phc_diffsmuw)
#                   elpd_diff se_diff
# ell_ph_diffsmuw   0.0       0.0   
# ell_phc_diffsmuw -0.8       1.0 

# Change in pH (DIW) and change in pH (CaCl2) equivalent predictors of Ellenberg
# R change. 

# Strength of relationship for pH:
# small weighted > large weighted > small unweighted > large unweighted
# for pH (CaCl2)
# large weighted ~= small weighted > small unweighted > large unweighted

# unweighted Ell R scores show very little relationship with pH change (95% overlaps 0)  



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
  left_join(BH_comb_nodupes) %>%
  left_join(select(CS_plot_atdep, REP_ID, Year, Ndep, Sdep) %>%
              mutate(Year = ifelse(Year == 2018, 2019, Year))) %>%
  left_join(cs_rainfall_diff) %>%
  left_join(cs_rainfall_averages) %>%
  left_join(PH_QA_diff) %>%
  left_join(rename(PH_long, Year1 = Year, Year1_pH = pH, 
                   Year1_pHCaCl2 = pH_CaCl2))
summary(ELL_pH)
mice::md.pattern(ELL_pH)

# run model with rain differences
ELL_pH_7807 <- filter(ELL_pH, Year != 2019) %>%
  select(-PHC, -PH_CACL2_SE, -Year1_pHCaCl2)
janitor::get_dupes(ELL_pH_7807, REP_ID, Time_period)
mice::md.pattern(ELL_pH_7807)


mod_data <- ELL_pH_7807 %>%
  mutate(YRnm = as.integer(as.factor(Year)),
         SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         Ell = WH_R_W) %>%
  filter(!BH %in% c(4,3,2,5,20,21,14,22,13,18,19,17)) %>%
  select(REP_ID, SQUARE, YRnm, Ell, Sdep,
         fieldseason_rain, Year1_pH, PH, PH_DIW_SE) %>%
  mutate(Sdep = as.numeric(scale(Sdep)), 
         fieldseason_rain = as.numeric(scale(fieldseason_rain)),
         ELL_SE = 0.16) %>%
  na.omit()
# 1085 obs

get_prior(bf(Ell | se(ELL_SE, sigma = TRUE) ~ PH*Year1_pH + (1|p|SQUARE) + ar(time = YRnm, gr = REP_ID),
             family = "student") +
            bf(PH | se(PH_DIW_SE, sigma = TRUE)  ~ Sdep + fieldseason_rain + (1|p|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), family = "student") + 
            set_rescor(FALSE),
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
                 data = mod_data, prior = mod_pr, cores = 6, iter = 5000)
summary(prior_mod)
pp_check(prior_mod, nsamples = 50, resp = "Ell")
pp_check(prior_mod, nsamples = 50, resp = "PH")
plot(conditional_effects(prior_mod))
# underestimates peak at 0 on average but gets close
stancode(prior_mod) # check if what I think is happening is actually happening

# model run
phell_mod <- brm(bf(Ell | se(ELL_SE, sigma = TRUE) ~ PH*Year1_pH + (1|p|SQUARE) + 
                      ar(time = YRnm, gr = REP_ID, cov = TRUE), family = "student") +
                   bf(PH | se(PH_DIW_SE, sigma = TRUE) ~ Sdep + fieldseason_rain + (1|p|SQUARE) + 
                        ar(time = YRnm, gr = REP_ID, cov = TRUE), family = "student") +
                   set_rescor(FALSE),  save_all_pars = TRUE,
                 data = mod_data, prior = mod_pr, cores = 6, iter = 5000,
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
                       data = mod_data, prior = mod_pr, cores = 6, iter = 5000,
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
                        data = mod_data, prior = mod_pr, cores = 6, iter = 5000,
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
                           data = mod_data, prior = mod_pr, cores = 6, iter = 5000,
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
                           data = mod_data, prior = mod_pr, cores = 6, iter = 5000,
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
                            data = mod_data, prior = mod_pr, cores = 6, iter = 5000,
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
                            data = mod_data, prior = mod_pr, cores = 6, iter = 5000,
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