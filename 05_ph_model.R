library(brms)
library(dplyr)
library(ggplot2)
library(tidyr)
theme_set(theme_classic())

str(PH_long)
str(LOI_long)
str(MOISTURE)
str(CS_plot_atdep)
str(cs_rainfall_stats)
str(BH_comb)

Soils_data <- full_join(PH_long, LOI_long) %>%
  left_join(MOISTURE) %>%
  left_join(mutate(CS_plot_atdep, 
                   Year = ifelse(Year == 2018, 2019, Year))) %>%
  left_join(cs_rainfall_stats) %>%
  left_join(BH_comb) %>%
  filter(!is.na(pH)) %>%
  select(REP_ID, Year, pH, BH_DESC, LOI, Ndep, Sdep,
         AVER_RAIN_8110, GWC = Moisture, rain_diff, pH_CaCl2) %>%
  mutate(SQUARE = sapply(strsplit(REP_ID, "[A-Z]"),"[",1))

mice::md.pattern(Soils_data)



ggplot(Soils_data, aes(x = LOI, y = pH)) + 
  geom_point(alpha = 0.2) +
  scale_x_log10()

ggplot(filter(Soils_data, !is.na(GWC)), 
       aes(x = LOI, y = pH, colour = GWC)) + 
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  facet_wrap(~BH_DESC)

# Broadleaved woodland only
Soils_wood <- Soils_data %>%
  filter(BH_DESC == "Broadleaved Mixed and Yew Woodland") %>%
  mutate(YR = as.factor(Year)) %>%
  mutate(YRnm = as.integer(YR)) %>%
  unique()

ggplot(Soils_wood, aes(x = LOI, y = pH)) + 
  geom_point(alpha = 0.5) +
  scale_x_log10()

ggplot(Soils_wood, aes(x = GWC, y = pH)) + 
  geom_point(alpha = 0.5)

ggplot(Soils_wood, aes(x = Sdep, y = pH)) + 
  geom_point(alpha = 0.5)

ggplot(Soils_wood, aes(x = AVER_RAIN_8110, y = pH)) + 
  geom_point(alpha = 0.5)

ggplot(Soils_wood, aes(x = rain_diff, y = pH)) + 
  geom_point(alpha = 0.5)

# model construction
get_prior(pH ~ LOI + GWC + Sdep + AVER_RAIN_8110 + rain_diff + (1|SQUARE) +
            ar(time = YRnm, gr = REP_ID),
          data = Soils_wood)
mod_pr <- c(prior(student_t(3,5,1), class = "Intercept"),
            prior(normal(0,0.5), class = "b"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"),
            prior(normal(0.5,0.1), class = "ar"))

Soils_wood_nona <- Soils_wood %>%
  select(REP_ID, SQUARE, Year, YRnm, 
         LOI, GWC, pH, Sdep, AVER_RAIN_8110) %>%
  na.omit()

ph_mod1 <- brm(pH ~ LOI + GWC + Sdep + AVER_RAIN_8110 + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_wood_nona, cores = 6, iter = 4000,
               save_all_pars = TRUE)
summary(ph_mod1)
plot(ph_mod1)
pp_check(ph_mod1)
ph_mod1 <- add_criterion(ph_mod1, "loo", moment_match = TRUE, reloo= TRUE)

ph_mod2 <- brm(pH ~ LOI + GWC + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_wood_nona, cores = 6, iter = 4000,
               save_all_pars = TRUE)
summary(ph_mod2)
plot(ph_mod2)
pp_check(ph_mod2)
ph_mod2 <- add_criterion(ph_mod2, "loo", moment_match = TRUE, reloo = TRUE)

loo_compare(ph_mod1, ph_mod2)

ph_mod3 <- brm(pH ~ LOI + GWC + AVER_RAIN_8110 + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_wood_nona, cores = 6, iter = 4000,
               save_all_pars = TRUE)
summary(ph_mod3)
plot(ph_mod3)
pp_check(ph_mod3)
ph_mod3 <- add_criterion(ph_mod3, "loo", moment_match = TRUE, reloo = TRUE)

ph_mod4 <- brm(pH ~ LOI + GWC + Sdep + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_wood_nona, cores = 6, iter = 4000,
               save_all_pars = TRUE)
summary(ph_mod4)
plot(ph_mod4)
pp_check(ph_mod4)
ph_mod4 <- add_criterion(ph_mod4, "loo", moment_match = TRUE, reloo = TRUE)

loo_compare(ph_mod1, ph_mod2, ph_mod3, ph_mod4)
loo_compare(ph_mod2, ph_mod3, ph_mod4)


ph_mod5 <- brm(pH ~ LOI + Sdep + AVER_RAIN_8110 + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_wood_nona, cores = 6, iter = 4000,
               save_all_pars = TRUE)
summary(ph_mod5)
plot(ph_mod5)
pp_check(ph_mod5)
ph_mod5 <- add_criterion(ph_mod5, "loo", moment_match = TRUE, reloo= TRUE)

loo_compare(ph_mod1, ph_mod2, ph_mod3, ph_mod4, ph_mod5)
#          elpd_diff se_diff
# ph_mod5  0.0       0.0   
# ph_mod1 -0.5       1.2   
# ph_mod4 -1.6       1.4   
# ph_mod3 -2.1       1.3   
# ph_mod2 -3.6       1.8   

# Full dataset
Soils_nona <- Soils_data %>%
  mutate(Wood = grepl("Woodland", BH_DESC),
         Bog = grepl("Bog", BH_DESC),
         YR = as.factor(Year)) %>%
  mutate(YRnm = as.integer(YR),
         Bog = replace_na(Bog, 0),
         Wood = replace_na(Wood, 0)) %>%
  select(REP_ID, SQUARE, Year, YRnm, Wood, Bog,
         LOI, GWC, pH, Sdep, AVER_RAIN_8110) %>%
  unique() %>%
  group_by(REP_ID, YRnm, SQUARE, Year, pH, LOI, GWC, Sdep, AVER_RAIN_8110) %>%
  summarise(Wood = max(Wood), Bog = max(Bog)) %>% ungroup() %>%
  na.omit() %>%
  mutate(across(LOI:AVER_RAIN_8110, scale))

ggplot(Soils_nona, aes(x = LOI, y = pH, colour = Bog)) + 
  geom_point(alpha = 0.5) +
  scale_x_log10()


ggplot(Soils_nona, aes(x = GWC, y = pH)) + 
  geom_point(alpha = 0.5)

ggplot(Soils_nona, aes(x = Sdep, y = pH, colour = Wood)) + 
  geom_point(alpha = 0.5)

ggplot(Soils_nona, aes(x = AVER_RAIN_8110, y = pH)) + 
  geom_point(alpha = 0.5)

ggplot(Soils_nona, aes(x = Wood, y = pH)) + 
  geom_jitter(alpha = 0.5, height = 0)

ggplot(Soils_nona, aes(x = Bog, y = pH)) + 
  geom_jitter(alpha = 0.5, height = 0)


# model construction
get_prior(pH ~ LOI*GWC + Sdep + AVER_RAIN_8110 + (1|SQUARE) +
            ar(time = YRnm, gr = REP_ID),
          data = Soils_nona)
mod_pr <- c(prior(student_t(3,5,1), class = "Intercept"),
            prior(normal(0,0.5), class = "b"),
            prior(student_t(3,0,1), class = "sd"),
            prior(student_t(3,0,1), class = "sigma"),
            prior(normal(0.5,0.1), class = "ar"))


ph_mod1 <- brm(pH ~ LOI*GWC + Sdep + AVER_RAIN_8110 + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_nona, cores = 4, iter = 5000,
               save_all_pars = TRUE, control = list(max_treedepth = 15),
               file = "Outputs/Models/ph_fulldataset_LOIGWC_Sdep_rain")
summary(ph_mod1)
plot(ph_mod1, ask = FALSE)
pp_check(ph_mod1)
ph_mod1 <- add_criterion(ph_mod1, "loo", moment_match = TRUE, reloo= TRUE)

ph_mod2 <- brm(pH ~ LOI*GWC + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_nona, cores = 4, iter = 5000,
               save_all_pars = TRUE,control = list(max_treedepth = 15),
               file = "Outputs/Models/ph_fulldataset_LOIGWC")
summary(ph_mod2)
plot(ph_mod2, ask = FALSE)
pp_check(ph_mod2)
ph_mod2 <- add_criterion(ph_mod2, "loo", moment_match = TRUE, reloo = TRUE)

loo_compare(ph_mod1, ph_mod2)

ph_mod3 <- brm(pH ~ LOI*GWC + AVER_RAIN_8110 + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_nona, cores = 4, iter = 5000,
               save_all_pars = TRUE,control = list(max_treedepth = 15),
               file = "Outputs/Models/ph_fulldataset_LOIGWC_rain")
summary(ph_mod3)
plot(ph_mod3, ask = FALSE)
pp_check(ph_mod3)
ph_mod3 <- add_criterion(ph_mod3, "loo", moment_match = TRUE, reloo = TRUE)

ph_mod4 <- brm(pH ~ LOI*GWC + Sdep + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_nona, cores = 4, iter = 5000,
               save_all_pars = TRUE,control = list(max_treedepth = 15),
               file = "Outputs/Models/ph_fulldataset_LOIGWC_Sdep")
summary(ph_mod4)
plot(ph_mod4, ask = FALSE)
pp_check(ph_mod4)
ph_mod4 <- add_criterion(ph_mod4, "loo", moment_match = TRUE, reloo = TRUE)

loo_compare(ph_mod1, ph_mod2, ph_mod3, ph_mod4)
loo_compare(ph_mod2, ph_mod3, ph_mod4)


ph_mod5 <- brm(pH ~ LOI + Sdep + AVER_RAIN_8110 + (1|SQUARE) +
                 ar(time = YRnm, gr = REP_ID), prior = mod_pr,
               data = Soils_nona, cores = 4, iter = 5000,
               save_all_pars = TRUE,control = list(max_treedepth = 15),
               file = "Outputs/Models/ph_fulldataset_LOI_Sdep_rain")
summary(ph_mod5)
plot(ph_mod5, ask = FALSE)
pp_check(ph_mod5)
ph_mod5 <- add_criterion(ph_mod5, "loo", moment_match = TRUE, reloo= TRUE)

loo_compare(ph_mod1, ph_mod2, ph_mod3, ph_mod4, ph_mod5)
