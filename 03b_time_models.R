# Summary statistics
library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_classic())
library(patchwork)
library(tidybayes)
library(emmeans)

# summary statistics
X_Ell %>% select(Year, contains("_R_")) %>%
  group_by(Year) %>%
  summarise(across(.fns = list(mean = ~mean(.x, na.rm = TRUE), 
                               sd = ~sd(.x, na.rm = TRUE), 
                               sem = ~ sd(.x, na.rm = TRUE)/sqrt(length(na.omit(.x))))))


X_Ell %>% select(Year, REP_ID, contains("_R_")) %>%
  pivot_longer(contains("_R_"),
               names_to = "Ell", values_to = "value") %>%
  pivot_wider(id_cols = REP_ID, names_from = c(Ell,Year), names_sep = "YR", 
              values_from = "value") %>%
  mutate(SM_R_W_diff9078 = SM_R_WYR1990 - SM_R_WYR1978,
         SM_R_W_diff9878 = SM_R_WYR1998 - SM_R_WYR1978,
         SM_R_W_diff9890 = SM_R_WYR1998 - SM_R_WYR1990,
         SM_R_W_diff0798 = SM_R_WYR2007 - SM_R_WYR1998,
         SM_R_W_diff1907 = SM_R_WYR2019 - SM_R_WYR2007,
         WH_R_W_diff9078 = WH_R_WYR1990 - WH_R_WYR1978,
         WH_R_W_diff9878 = WH_R_WYR1998 - WH_R_WYR1978,
         WH_R_W_diff9890 = WH_R_WYR1998 - WH_R_WYR1990,
         WH_R_W_diff0798 = WH_R_WYR2007 - WH_R_WYR1998,
         WH_R_W_diff1907 = WH_R_WYR2019 - WH_R_WYR2007,
         SM_R_UW_diff9078 = SM_R_UWYR1990 - SM_R_UWYR1978,
         SM_R_UW_diff9878 = SM_R_UWYR1998 - SM_R_UWYR1978,
         SM_R_UW_diff9890 = SM_R_UWYR1998 - SM_R_UWYR1990,
         SM_R_UW_diff0798 = SM_R_UWYR2007 - SM_R_UWYR1998,
         SM_R_UW_diff1907 = SM_R_UWYR2019 - SM_R_UWYR2007,
         WH_R_UW_diff9078 = WH_R_UWYR1990 - WH_R_UWYR1978,
         WH_R_UW_diff9878 = WH_R_UWYR1998 - WH_R_UWYR1978,
         WH_R_UW_diff9890 = WH_R_UWYR1998 - WH_R_UWYR1990,
         WH_R_UW_diff0798 = WH_R_UWYR2007 - WH_R_UWYR1998,
         WH_R_UW_diff1907 = WH_R_UWYR2019 - WH_R_UWYR2007
  ) %>% select(contains("diff")) %>%
  summarise(across(.fns = list(mean = ~mean(.x, na.rm = TRUE), 
                               sem = ~ sd(.x, na.rm = TRUE)/sqrt(length(na.omit(.x)))))) %>%
  pivot_longer(everything(),names_to = c("Variable", "Year"), names_sep = "_diff") %>%
  print(n = 40)

summary(PH)
PH %>% select(-REP_ID) %>%
  summarise(across(.fns = list(mean = ~mean(.x, na.rm = TRUE), 
                               sem = ~ sd(.x, na.rm = TRUE)/sqrt(length(na.omit(.x))))))
PH_C_0719 <- na.omit(PH$PHC_2019 - PH$PHC_2007)
mean(PH_C_0719)
sd(PH_C_0719)/sqrt(length(PH_C_0719))


# Bayesian models
mod_data <- mutate(X_Ell, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = SM_R_W,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR)) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL) #%>%
# na.omit()



# Year models 
get_prior(ELL ~ YR - 1 + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID), data = mod_data)
mod_pr <- c(prior(normal(5,0.5), class = "b"),
            prior(normal(0.2,0.1), class = "ar"),
            prior(student_t(3, 0, 1), class = "sd"))

# weighted 4m2 Ellenberg R
sm_r_w_mod <- brm(ELL ~ YR - 1 + (1|SERIES_NUM) +
                    ar(time = YRnm, gr = REP_ID), prior = mod_pr,
                  data = mod_data, cores = 8, iter = 8000,
                  file = paste0("Outputs/Models/Year_models/SM_R_W-brms-", Sys.Date()))

summary(sm_r_w_mod)
pp_check(sm_r_w_mod)
plot(conditional_effects(sm_r_w_mod))

# weighted 200m2 Ellenberg R
mod_data <- mutate(X_Ell, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = WH_R_W,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR)) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL) %>%
  na.omit()
wh_r_w_mod <- update(sm_r_w_mod, newdata = mod_data, cores = 8, iter = 8000,
                     file = paste0("Outputs/Models/Year_models/WH_R_W-brms-", Sys.Date()))
summary(wh_r_w_mod)
plot(wh_r_w_mod)
pp_check(wh_r_w_mod)
pp_check(wh_r_w_mod, type = "stat_2d")
pp_check(wh_r_w_mod, type = "stat_grouped", group = "YR")
pr <- predict(wh_r_w_mod, summary = FALSE)
str(pr)
bayesplot::ppc_stat_grouped(as.vector(na.omit(mod_data$ELL)), pr,
                            group = as.character(
                              as.data.frame(mod_data)[!is.na(mod_data$ELL),"YR"]))
plot(conditional_effects(wh_r_w_mod))

# unweighted 200m2 Ellenberg R
mod_data <- mutate(mod_data, ELL = WH_R_UW)
wh_r_uw_mod <- update(sm_r_w_mod, newdata = mod_data, cores = 8, iter = 8000,
                      file = paste0("Outputs/Models/Year_models/WH_R_UW-brms-", Sys.Date()))
summary(wh_r_uw_mod)
plot(wh_r_uw_mod)
pp_check(wh_r_uw_mod)
pp_check(wh_r_uw_mod, type = "stat_2d")
pp_check(wh_r_uw_mod, type = "stat_grouped", group = "YR")
plot(conditional_effects(wh_r_uw_mod))

# unweighted 4m2 Ellenberg R
mod_data <- mutate(mod_data, ELL = SM_R_UW)
sm_r_uw_mod <- update(sm_r_w_mod, newdata = mod_data, cores = 8, iter = 8000,
                      file = paste0("Outputs/Models/Year_models/SM_R_UW-brms-", Sys.Date()))
summary(sm_r_uw_mod)
plot(sm_r_uw_mod)
pp_check(sm_r_uw_mod)
pp_check(sm_r_uw_mod, type = "stat_2d")
pp_check(sm_r_uw_mod, type = "stat_grouped", group = "YR")
plot(conditional_effects(sm_r_uw_mod))


wh_w <- wh_r_w_mod %>%
  emmeans( ~ YR) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value)) +
  stat_lineribbon(alpha = 1/4, fill = "#2F7ECE") +
  scale_y_continuous(limits = c(4.7,5.6))+
  labs(y = bquote("Cover weighted Ellenberg R 200m"^2))
wh_uw <- wh_r_uw_mod %>%
  emmeans( ~ YR) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value)) +
  stat_lineribbon(alpha = 1/4, fill = "#2F7ECE") +
  scale_y_continuous(limits = c(4.7,5.6)) +
  labs(y = bquote("Unweighted Ellenberg R 200m"^2))
sm_w <- sm_r_w_mod %>%
  emmeans( ~ YR) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value)) +
  stat_lineribbon(alpha = 1/4, fill = "#2F7ECE") +
  scale_y_continuous(limits = c(4.7,5.6)) +
  labs(y = bquote("Cover weighted Ellenberg R 4m"^2))
sm_uw <- sm_r_uw_mod %>%
  emmeans( ~ YR) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value)) +
  stat_lineribbon(alpha = 1/4, fill = "#2F7ECE") +
  scale_y_continuous(limits = c(4.7,5.6)) +
  labs(y = bquote("Unweighted Ellenberg R 4m"^2))

wh_w + wh_uw + sm_w + sm_uw + plot_layout(ncol = 2)
ggsave("Ellenberg R versus time model results.png",
       path = "Outputs/Models/Year_models/", width = 20, height = 20, units = "cm")

# pH models ####
summary(PH_long)
mod_data <- PH_long %>%
  mutate(Year = ifelse(Year == 2016, 2019, Year)) %>%
  mutate(SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         PH = pH,
         YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR)) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, PH) %>%
  na.omit()

get_prior(PH ~ YR - 1 + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID), data = mod_data)
mod_pr <- c(prior(normal(5,0.5), class = "b"),
            prior(normal(0.2,0.1), class = "ar"),
            prior(student_t(3, 0, 2.5), class = "sd"))

# pH in water 
ph_mod <- brm(PH ~ YR - 1 + (1|SERIES_NUM) +
                ar(time = YRnm, gr = REP_ID),
              data = mod_data, cores = 8, iter = 8000,
              file = paste0("Outputs/Models/Year_models/ph_mod-brms-", Sys.Date()))
summary(ph_mod)
plot(ph_mod)
pp_check(ph_mod)
pp_check(ph_mod, type = "stat_2d")
pp_check(ph_mod, type = "stat_grouped", group = "YR")
pr <- predict(ph_mod, summary = FALSE)
str(pr)
bayesplot::ppc_stat_grouped(as.vector(na.omit(mod_data$PH)), pr,
                            group = as.character(
                              as.data.frame(mod_data)[!is.na(mod_data$PH),"YR"]))
plot(conditional_effects(ph_mod))

ph_pl <- ph_mod %>%
  emmeans( ~ YR) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value)) +
  stat_lineribbon(alpha = 1/4, fill = "#2F7ECE") +
  scale_y_continuous(limits = c(5.1,6.1)) +
  labs(y = "pH")


# pH in CaCl2
mod_data <-PH_long %>%
  mutate(Year = ifelse(Year == 2016, 2019, Year)) %>%
  mutate(SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
         PH = pH_CaCl2,
         YR = as.factor(Year)) %>%
  ungroup() %>%
  select(REP_ID, SERIES_NUM, YR, PH) %>%
  na.omit() %>% droplevels()%>%
  mutate(YRnm = as.integer(YR))

phc_mod <- update(ph_mod, newdata = mod_data,
                  cores = 8, iter = 8000,
                  file = paste0("Outputs/Models/Year_models/ph_c_mod-brms",Sys.Date()))
summary(phc_mod)
plot(phc_mod)
pp_check(phc_mod)
pp_check(phc_mod, type = "stat_2d")
pp_check(phc_mod, type = "stat_grouped", group = "YR")
pr <- predict(phc_mod, summary = FALSE)
str(pr)
bayesplot::ppc_stat_grouped(as.vector(na.omit(mod_data$PH)), pr,
                            group = as.character(
                              as.data.frame(mod_data)[!is.na(mod_data$PH),"YR"]))
plot(conditional_effects(phc_mod))

ph_c_pl <- phc_mod %>%
  emmeans( ~ YR) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value)) +
  stat_lineribbon(alpha = 1/4, fill = "#2F7ECE") +
  scale_y_continuous(limits = c(4.5,5.6)) +
  scale_x_continuous(limits = c(1978,2019)) +
  labs(y = bquote("pH (CaCl"[2]*")"))

ph_pl + ph_c_pl
ggsave("pH by year model results.png",
       path = "Outputs/Models/Year_models/", width =20, height = 10, units = "cm")

wh_w + wh_uw + sm_w + sm_uw + 
  ph_pl + ph_c_pl + plot_layout(ncol = 2)
ggsave("Ellenberg R and pH by year model results.png",
       path = "Outputs/Models/Year_models/", width =20, height = 30, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# Habitat ####
# combined data
Comb <- full_join(X_Ell, PH_long) %>% left_join(BH_IMP) %>%
  filter(!is.na(Management))
str(Comb)
mice::md.pattern(select(Comb, Year, REP_ID, Management, contains("_R_"),
                        pH, pH_CaCl2))
janitor::get_dupes(Comb, REP_ID, Year) # no dupes

Comb %>%
  select(Year, Management, contains("_R_"),
         pH, pH_CaCl2) %>%
  group_by(Year, Management) %>%
  summarise(across(.fns = list(LLmean = ~ mean(.x, na.rm = TRUE),
                               LLsem = ~ sd(.x, na.rm = TRUE)/sqrt(length(na.omit(.x)))))) %>%
  pivot_longer(cols = SM_R_W_LLmean:pH_CaCl2_LLsem,
               names_to = c("Variable","Type"),
               names_sep = "LL") %>%
  pivot_wider(names_from = Year, values_from = value) %>%
  arrange(Variable)

mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = SM_R_W,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL, Management) %>%
  na.omit()

# Year models 
get_prior(ELL ~ YR*Management + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID), data = mod_data)
mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(5,0.5), class = "Intercept"),
            prior(normal(0.5,0.2), class = "ar"),
            prior(student_t(3, 0, 1), class = "sd"))

# weighted 4m2 Ellenberg R
sm_r_w_hab_mod <- brm(ELL ~ YR*Management + (1|SERIES_NUM) +
                        ar(time = YRnm, gr = REP_ID), prior = mod_pr,
                      data = mod_data, cores = 4, iter = 4000,
                      file = paste0("Outputs/Models/Year_models/SM_R_W_HAB-brms-", Sys.Date()))
summary(sm_r_w_hab_mod)
plot(sm_r_w_hab_mod)
pp_check(sm_r_w_hab_mod)
pp_check(sm_r_w_hab_mod, type = "stat_2d")
# pr <- predict(sm_r_w_hab_mod, summary = FALSE)
# str(pr)
# bayesplot::ppc_stat_grouped(as.vector(na.omit(mod_data$ELL)), pr,
#                             group = as.character(
#                               as.data.frame(mod_data)[!is.na(mod_data$ELL),"YR"]))
plot(conditional_effects(sm_r_w_hab_mod), ask = FALSE)


# unweighted 4m2 plot
mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = SM_R_UW,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL, Management) %>%
  na.omit()

sm_r_uw_hab_mod <- update(sm_r_w_hab_mod, newdata = mod_data, cores = 4, iter = 4000,
                          file = paste0("Outputs/Models/Year_models/SM_R_UW_HAB-brms",Sys.Date()))
summary(sm_r_uw_hab_mod)
plot(sm_r_uw_hab_mod)
pp_check(sm_r_uw_hab_mod)
pp_check(sm_r_uw_hab_mod, type = "stat_2d")
# pr <- predict(sm_r_uw_hab_mod, summary = FALSE)
# str(pr)
# bayesplot::ppc_stat_grouped(as.vector(na.omit(mod_data$ELL)), pr,
#                             group = as.character(
#                               as.data.frame(mod_data)[!is.na(mod_data$ELL),"YR"]))
plot(conditional_effects(sm_r_uw_hab_mod), ask = FALSE)


# unweighted 200m2 plot
mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = WH_R_UW,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL, Management) %>%
  na.omit()

wh_r_uw_hab_mod <- update(sm_r_w_hab_mod, newdata = mod_data, cores = 4, iter = 4000,
                          file = paste0("Outputs/Models/Year_models/WH_R_UW_HAB-brms",Sys.Date()))
summary(wh_r_uw_hab_mod)
plot(wh_r_uw_hab_mod)
pp_check(wh_r_uw_hab_mod)
pp_check(wh_r_uw_hab_mod, type = "stat_2d")
# pr <- predict(wh_r_uw_hab_mod, summary = FALSE)
# str(pr)
# bayesplot::ppc_stat_grouped(as.vector(na.omit(mod_data$ELL)), pr,
#                             group = as.character(
#                               as.data.frame(mod_data)[!is.na(mod_data$ELL),"YR"]))
plot(conditional_effects(wh_r_uw_hab_mod), ask = FALSE)

# unweighted 200m2 plot
mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = WH_R_W,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL, Management) %>%
  na.omit()

wh_r_w_hab_mod <- update(sm_r_w_hab_mod, newdata = mod_data, cores = 4, iter = 4000,
                         file = paste0("Outputs/Models/Year_models/WH_R_W_HAB-brms",Sys.Date()))
summary(wh_r_w_hab_mod)
plot(wh_r_w_hab_mod)
pp_check(wh_r_w_hab_mod)
pp_check(wh_r_w_hab_mod, type = "stat_2d")
# pr <- predict(wh_r_w_hab_mod, summary = FALSE)
# str(pr)
# bayesplot::ppc_stat_grouped(as.vector(na.omit(mod_data$ELL)), pr,
#                             group = as.character(
#                               as.data.frame(mod_data)[!is.na(mod_data$ELL),"YR"]))
plot(conditional_effects(wh_r_w_hab_mod), ask = FALSE)


# nice conditional effect plots
wh_w <- wh_r_w_hab_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value, 
             fill = Management, colour = Management)) +
  stat_lineribbon(alpha = 1/4) +
  scale_y_continuous(limits = c(4,6))+
  labs(y = bquote("Cover weighted Ellenberg R 200m"^2))
wh_uw <- wh_r_uw_hab_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value, 
             fill = Management, colour = Management)) +
  stat_lineribbon(alpha = 1/4) +
  scale_y_continuous(limits = c(4,6)) +
  labs(y = bquote("Unweighted Ellenberg R 200m"^2))
sm_w <- sm_r_w_hab_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value,
             fill = Management, colour = Management)) +
  stat_lineribbon(alpha = 1/4) +
  scale_y_continuous(limits = c(4,6)) +
  labs(y = bquote("Cover weighted Ellenberg R 4m"^2))
sm_uw <- sm_r_uw_hab_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value,
             colour = Management, fill = Management)) +
  stat_lineribbon(alpha = 1/4) +
  scale_y_continuous(limits = c(4,6)) +
  labs(y = bquote("Unweighted Ellenberg R 4m"^2))

wh_w + wh_uw + sm_w + sm_uw + plot_layout(guides = "collect") 
ggsave("Ellenberg R versus time by Management variable y limits.png",
       path = "Outputs/Models/Year_models/",
       width = 20, height = 20, units = "cm")

# Ell w/meas error ####

ELL_QA_year <- data.frame(
  Year = c(1978,1990,1998,2007,2019),
  SM_UW_SE_NORM = c(0.240,0.240,0.229,0.308,0.320),
  SM_W_SE_NORM = c(0.396,0.396,0.312,0.468,0.480),
  WH_UW_SE_NORM = c(0.204,0.204,0.198,0.236,0.181),
  WH_W_SE_NORM = c(0.298,0.298,0.239,0.319,0.419)
)

Comb <- left_join(Comb, ELL_QA_year)


mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = SM_R_W,
                   ELL_SE = SM_W_SE_NORM,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL, ELL_SE, Management) %>%
  na.omit()

# Year models 
get_prior(ELL | mi(ELL_SE) ~ YR*Management + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID), data = mod_data)
mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(5,0.5), class = "Intercept"),
            prior(normal(0.5,0.2), class = "ar"),
            prior(student_t(3, 0, 1), class = "sd"))

# weighted 4m2 Ellenberg R
sm_r_w_hab_me_mod <- brm(ELL | mi(ELL_SE) ~ YR*Management + (1|SERIES_NUM) +
                           ar(time = YRnm, gr = REP_ID), prior = mod_pr,
                         data = mod_data, cores = 4, iter = 6000,
                         file = paste0("Outputs/Models/Year_models/SM_R_W_HAB_me-brms-", Sys.Date()))
summary(sm_r_w_hab_me_mod)
plot(sm_r_w_hab_me_mod, ask = FALSE)
pp_check(sm_r_w_hab_me_mod)
pp_check(sm_r_w_hab_me_mod, type = "stat_2d")
plot(conditional_effects(sm_r_w_hab_me_mod), ask = FALSE)


# unweighted 4m2 plot
mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = SM_R_UW,
                   ELL_SE = SM_UW_SE_NORM,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL, ELL_SE, Management) %>%
  na.omit()

sm_r_uw_hab_me_mod <- update(sm_r_w_hab_me_mod, newdata = mod_data, 
                             cores = 4, iter = 6000,
                             file = paste0("Outputs/Models/Year_models/SM_R_UW_HAB_me-brms",Sys.Date()))
summary(sm_r_uw_hab_me_mod)
plot(sm_r_uw_hab_me_mod, ask = FALSE)
pp_check(sm_r_uw_hab_me_mod)
pp_check(sm_r_uw_hab_me_mod, type = "stat_2d")
plot(conditional_effects(sm_r_uw_hab_me_mod), ask = FALSE)


# unweighted 200m2 plot
mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = WH_R_UW,
                   ELL_SE = WH_UW_SE_NORM,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL, ELL_SE, Management) %>%
  na.omit()

wh_r_uw_hab_me_mod <- update(sm_r_w_hab_me_mod, newdata = mod_data, 
                             cores = 4, iter = 6000,
                             file = paste0("Outputs/Models/Year_models/WH_R_UW_HAB_me-brms",Sys.Date()))
summary(wh_r_uw_hab_me_mod)
plot(wh_r_uw_hab_me_mod, ask = FALSE)
pp_check(wh_r_uw_hab_me_mod)
pp_check(wh_r_uw_hab_me_mod, type = "stat_2d")
plot(conditional_effects(wh_r_uw_hab_me_mod), ask = FALSE)

# unweighted 200m2 plot
mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   ELL = WH_R_W,
                   ELL_SE = WH_R_SE_NORM,
                   YR = as.factor(Year)) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, ELL, ELL_SE, Management) %>%
  na.omit()

wh_r_w_hab_me_mod <- update(sm_r_w_hab_me_mod, newdata = mod_data, 
                            cores = 4, iter = 6000,
                            file = paste0("Outputs/Models/Year_models/WH_R_W_HAB_me-brms",Sys.Date()))
summary(wh_r_w_hab_me_mod)
plot(wh_r_w_hab_me_mod, ask = FALSE)
pp_check(wh_r_w_hab_me_mod)
pp_check(wh_r_w_hab_me_mod, type = "stat_2d")
plot(conditional_effects(wh_r_w_hab_me_mod), ask = FALSE)


# nice conditional effect plots
wh_w <- wh_r_w_hab_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value, 
             fill = Management, colour = Management)) +
  stat_lineribbon(alpha = 1/4) +
  scale_y_continuous(limits = c(4,6))+
  labs(y = bquote("Cover weighted Ellenberg R 200m"^2))
wh_uw <- wh_r_uw_hab_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value, 
             fill = Management, colour = Management)) +
  stat_lineribbon(alpha = 1/4) +
  scale_y_continuous(limits = c(4,6)) +
  labs(y = bquote("Unweighted Ellenberg R 200m"^2))
sm_w <- sm_r_w_hab_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value,
             fill = Management, colour = Management)) +
  stat_lineribbon(alpha = 1/4) +
  scale_y_continuous(limits = c(4,6)) +
  labs(y = bquote("Cover weighted Ellenberg R 4m"^2))
sm_uw <- sm_r_uw_hab_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  mutate(Year = as.numeric(as.character(YR))) %>%
  ggplot(aes(x = Year, y = .value,
             colour = Management, fill = Management)) +
  stat_lineribbon(alpha = 1/4) +
  scale_y_continuous(limits = c(4,6)) +
  labs(y = bquote("Unweighted Ellenberg R 4m"^2))

wh_w + wh_uw + sm_w + sm_uw + plot_layout(guides = "collect")
ggsave("Ellenberg R versus time by Management with measurement error variable y limits.png",
       path = "Temp/", width = 30, height = 30, units = "cm")

# pH models ####
mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   PH = pH,
                   YR = as.factor(Year)) %>%
  filter(Year != 2016) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, PH, Management) %>%
  na.omit() %>% droplevels()

# Year models 
get_prior(PH ~ YR*Management + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID), data = mod_data)
mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(5,0.5), class = "Intercept"),
            prior(normal(0.5,0.2), class = "ar"),
            prior(student_t(3, 0, 1), class = "sd"))

#  pH in water
ph_hab_mod <- brm(PH ~ YR*Management + (1|SERIES_NUM) +
                    ar(time = YRnm, gr = REP_ID), prior = mod_pr,
                  data = mod_data, cores = 4, iter = 4000,
                  file = paste0("Outputs/Models/Year_models/PH_hab-brms-", Sys.Date()))
summary(ph_hab_mod)
plot(ph_hab_mod)
pp_check(ph_hab_mod)
pp_check(ph_hab_mod, type = "stat_2d")
pr <- predict(ph_hab_mod, summary = FALSE)
str(pr)
bayesplot::ppc_stat_grouped(as.vector(na.omit(mod_data$PH)), pr,
                            group = as.character(
                              as.data.frame(mod_data)[!is.na(mod_data$PH),"YR"]))
plot(conditional_effects(ph_hab_mod), ask = FALSE)

ph_pl <- ph_hab_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  filter(Management != "Other") %>%
  mutate(Year = as.numeric(as.character(YR)),
         Management = recode(Management, "Management" = "Other")) %>%
  ggplot(aes(x = Year, y = .value)) +
  stat_lineribbon(alpha = 1/4, fill = "#2F7ECE") +
  facet_wrap(~Management) +
  labs(y = "pH")
ph_pl
ggsave("pH versus time by Management.png", path = "Outputs/Models/Year_models/",
       width = 20, height = 18, units = "cm")

# pH in CaCl2
mod_data <- mutate(Comb, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   PH = pH_CaCl2,
                   YR = as.factor(Year)) %>%
  filter(Year != 2016) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, PH, Management) %>%
  na.omit() %>% droplevels()

ph_c_hab_mod <- update(ph_hab_mod, newdata = mod_data, cores = 4, iter = 4000,
                       file = paste0("Outputs/Models/Year_models/PH_C_HAB-brms",Sys.Date()))
summary(ph_c_hab_mod)
plot(ph_c_hab_mod)
pp_check(ph_c_hab_mod)
pp_check(ph_c_hab_mod, type = "stat_2d")
pr <- predict(ph_c_hab_mod, summary = FALSE)
str(pr)
bayesplot::ppc_stat_grouped(as.vector(na.omit(mod_data$PH)), pr,
                            group = as.character(
                              as.data.frame(mod_data)[!is.na(mod_data$PH),"YR"]))
plot(conditional_effects(ph_c_hab_mod), ask = FALSE)

# pH w/meas error ####
PH_QA <- data.frame(Year = c(2007,2019),
                    PH_SE = c(0.252,0.165),
                    PHC_SE = c(0.226,0.119))

Comb_PH_QA <- left_join(Comb, PH_QA) %>%
  filter(Year > 2005)

mod_data <- mutate(Comb_PH_QA, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   PH = pH,
                   YR = as.factor(Year)) %>%
  filter(Year != 2016) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, PH, PH_SE, Management) %>%
  na.omit() %>% droplevels()

# Year models 
get_prior(PH | mi(PH_SE) ~ YR*Management + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID), data = mod_data)
mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(5,0.5), class = "Intercept"),
            prior(normal(0.5,0.2), class = "ar"),
            prior(student_t(3, 0, 1), class = "sd"))

#  pH in water
ph_hab_me_mod <- brm(PH | mi(PH_SE) ~ YR*Management + (1|SERIES_NUM) +
                       ar(time = YRnm, gr = REP_ID), prior = mod_pr,
                     data = mod_data, cores = 4, iter = 4000,
                     file = paste0("Outputs/Models/Year_models/PH_hab_me-brms-", Sys.Date()))
summary(ph_hab_me_mod)
plot(ph_hab_me_mod)
pp_check(ph_hab_me_mod)
pp_check(ph_hab_me_mod, type = "stat_2d")
plot(conditional_effects(ph_hab_me_mod), ask = FALSE)

ph_pl <- ph_hab_me_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  filter(Management != "Other") %>%
  mutate(Year = as.numeric(as.character(YR)),
         Management = recode(Management, "Management" = "Other")) %>%
  ggplot(aes(x = Year, y = .value)) +
  stat_lineribbon(alpha = 1/4, fill = "#2F7ECE") +
  facet_wrap(~Management) +
  labs(y = "pH (DIW)")

# pH in CaCl2
mod_data <- mutate(Comb_PH_QA, 
                   SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1),
                   PH = pH_CaCl2,
                   PH_SE = PHC_SE,
                   YR = as.factor(Year)) %>%
  filter(Year != 2016) %>%
  ungroup() %>%
  mutate(YRnm = as.integer(YR),
         Management = replace_na(Management, "Other")) %>%
  select(REP_ID, SERIES_NUM, YR, YRnm, PH, Management) %>%
  na.omit() %>% droplevels()

ph_c_hab_me_mod <- update(ph_hab_me_mod, newdata = mod_data, cores = 4, iter = 4000,
                          file = paste0("Outputs/Models/Year_models/PH_C_HAB_me-brms",Sys.Date()))
summary(ph_c_hab_me_mod)
plot(ph_c_hab_me_mod)
pp_check(ph_c_hab_me_mod)
pp_check(ph_c_hab_me_mod, type = "stat_2d")
plot(conditional_effects(ph_c_hab_me_mod), ask = FALSE)

phc_pl <- ph_c_hab_me_mod %>%
  emmeans( ~ YR | Management) %>%
  gather_emmeans_draws() %>%
  filter(Management != "Other") %>%
  mutate(Year = as.numeric(as.character(YR)),
         Management = recode(Management, "Management" = "Other")) %>%
  ggplot(aes(x = Year, y = .value)) +
  stat_lineribbon(alpha = 1/4, fill = "#2F7ECE") +
  facet_wrap(~Management) +
  labs(y = bquote("pH CaCl"[2]))
ph_pl + phc_pl
ggsave("pH against year with measurement error.png",
       path = "Outputs/Models/Year_models/",
       width = 20, height = 18, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####
# pH against rainfall ####
str(Comb)
Comb_rain <- cs_survey_rainfall %>%
  mutate(Year = as.numeric(Year)) %>%
  full_join(Comb) %>%
  full_join(cs_loc_rain30y) %>%
  full_join(cs_loc_sumrain30y) %>%
  mutate(rain_diff_sum = mean_rainfall - SUM/3)

coplot(pH ~ mean_rainfall|RAIN_8110, data = Comb_rain)

ggplot(Comb_rain, aes(x = mean_rainfall, y = pH, colour = Year)) +
  geom_point()

ggplot(Comb_rain, aes(x = RAIN_8110, y = pH, colour = Year)) +
  geom_point()

Comb_rain %>%
  ggplot(aes(x=rain_diff_sum, y = pH, colour = YEAR)) +
  geom_point()


# Bayesian model
mod_data <- Comb_rain %>% ungroup() %>%
  select(Year, REP_ID, mean_rainfall, pH, Management, rain_diff_sum, RAIN_8110) %>%
  mutate(YR = as.factor(Year),
         SERIES_NUM = sapply(strsplit(REP_ID, "[A-Z]"),"[",1)) %>%
  mutate(YRnm = as.integer(YR)) %>%
  na.omit()

get_prior(pH ~ Year + Management + rain_diff_sum + RAIN_8110 + (1|SERIES_NUM) +
            ar(time = YRnm, gr = REP_ID), data = mod_data)
mod_pr <- c(prior(normal(0,0.5), class = "b"),
            prior(normal(5,0.5), class = "Intercept"),
            prior(normal(0.5,0.2), class = "ar"),
            prior(student_t(3, 0, 1), class = "sd"))
ph_rain_mod <- brm(pH ~ YR + Management + rain_diff_sum*RAIN_8110 + (1|SERIES_NUM) +
                     ar(time = YRnm, gr = REP_ID), data = mod_data,
                   prior = mod_pr, cores = 4, iter = 4000,
                   file = paste0("Outputs/Models/Year_models/ph_rain_mod_",Sys.Date()))
# took forever to run and failed miserably

ph_rain_mod_nohab <- brm(pH ~ YR + rain_diff_sum*RAIN_8110 + (1|SERIES_NUM) +
                           ar(time = YRnm, gr = REP_ID), data = mod_data,
                         prior = mod_pr, cores = 4, iter = 4000,
                         file = paste0("Outputs/Models/Year_models/ph_rain_mod_",Sys.Date()))
