# Summary statistics
library(brms)

mod_data <- mutate(X_Ell, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "X"),"[",1),
                   ELL = SM_R_W,
                   YR = as.factor(Year)) %>%
  mutate(YRnm = as.integer(YR))

# Year models 
get_prior(ELL ~ YR - 1 + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID), data = mod_data)
mod_pr <- c(prior(normal(5,0.5), class = "b"),
            prior(normal(0.2,0.1), class = "ar"),
            prior(student_t(3, 0, 1), class = "sd"))

# weighted 4m2 Ellenberg R
sm_r_w_mod <- brm(ELL ~ YR - 1 + (1|SERIES_NUM) +
                    ar(time = YRnm, gr = REP_ID),
                  data = mod_data, cores = 8, iter = 8000,
                  file = paste0("Outputs/Models/Year_models/SM_R_W-brms-", Sys.Date()))

summary(sm_r_w_mod)
pp_check(sm_r_w_mod)
# plot(conditional_effects(sm_r_w_mod))

# weighted 200m2 Ellenberg R
mod_data <- mutate(mod_data, ELL = WH_R_W)
wh_r_w_mod <- update(sm_r_w_mod, newdata = mod_data, cores = 8, iter = 8000,
                     file = paste0("~/Temp/WH_R_W-brms-", Sys.Date()))
summary(wh_r_w_mod)
pp_check(wh_r_w_mod)
# plot(conditional_effects(wh_r_w_mod))

# unweighted 200m2 Ellenberg R
mod_data <- mutate(mod_data, ELL = WH_R_UW)
wh_r_uw_mod <- update(sm_r_uw_mod, newdata = mod_data, cores = 8, iter = 8000,
                      file = paste0("~/Temp/WH_R_UW-brms-", Sys.Date()))
summary(wh_r_uw_mod)
pp_check(wh_r_uw_mod)
# plot(conditional_effects(wh_r_uw_mod))

# unweighted 4m2 Ellenberg R
mod_data <- mutate(mod_data, ELL = SM_R_UW)
sm_r_uw_mod <- update(sm_r_uw_mod, newdata = mod_data, cores = 8, iter = 8000,
                      file = paste0("~/Temp/SM_R_UW-brms-", Sys.Date()))
summary(sm_r_uw_mod)
pp_check(sm_r_uw_mod)
# plot(conditional_effects(sm_r_uw_mod))