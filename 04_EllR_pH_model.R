# Ellenberg R and pH model


library(brms)
library(ggplot2)
theme_set(theme_classic())

# taking ph_ell_wide from 02 script
ph_ell_long$SQNUM <- sapply(strsplit(ph_ell_long$REP_ID, "X"), "[",1)
ph_ell_long <- filter(ph_ell_long, year != 2018) %>%
  mutate(YR = as.factor(year))

# Non-linear model linking pH and Ellenberg R
# Assymetrical sigmoidal curve, with impact of pH varying by year
# c1 is the max y (minus c5)
# c2 is the rate of increase for the curve
# c3 is the x at curve midpoint
# c4 is the amount of asymmetry
# c5 is the min y
# c6 is the random effect(square shifts whole curve up/down)

get_prior(bf(Ell_R ~ c1/((1 + exp(-c2*(Soil_pH - c3)))^c4) + c5, 
             c5 ~ (1|SQNUM), c1 ~ 1, c2 + c3 + c4 ~ YR, nl= TRUE),
          data = ph_ell_long)

pr <- c(prior(normal(4,0.5), nlpar = "c1"),
        prior(normal(1.5,0.5), nlpar = "c2"),
        prior(normal(0,0.1), nlpar = "c2", coef = "YR1998"),
        prior(normal(0,0.1), nlpar = "c2", coef = "YR2007"),
        prior(normal(3.5,1), nlpar = "c3", coef = "Intercept"),
        prior(normal(0,0.1), nlpar = "c3", coef = "YR1998"),
        prior(normal(0,0.1), nlpar = "c3", coef = "YR2007"),
        prior(gamma(2,2), nlpar = "c4", lb = 0),
        prior(normal(0,0.1), nlpar = "c4", coef = "YR1998"),
        prior(normal(0,0.1), nlpar = "c4", coef = "YR2007"),
        prior(normal(2.5,0.5), nlpar = "c5", coef = "Intercept"),
        prior(student_t(3,0,.1), nlpar = "c5", group = "SQNUM", class = "sd"))

# prior simulation
mod_pr_only <- brm(bf(Ell_R ~ c1/((1 + exp(-c2*(Soil_pH - c3)))^c4) + c5, 
                      c5 ~ (1|SQNUM), c1 ~ 1, c2 + c3 + c4 ~ YR, nl= TRUE),
                   data = ph_ell_long, prior = pr, cores = 2, chains = 2,
                   sample_prior = "only")
summary(mod_pr_only)
plot(mod_pr_only)

pp_check(mod_pr_only)

pred_modpr <- predict(mod_pr_only)

plot(na.omit(ph_ell_long)$Ell_R, pred_modpr[,1])
