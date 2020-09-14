# Ellenberg R and pH model
library(brms)
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
