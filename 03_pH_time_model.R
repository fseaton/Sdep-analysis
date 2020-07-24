library(nlme)
library(emmeans)
library(dplyr)
library(car)
library(tidyr)

# data manipulation
CS78_PH$REP_ID <- paste(CS78_PH$SQUARE_NUM,CS78_PH$REP_NUM, sep = "X")
CS98_PH$REP_ID <- paste(CS98_PH$SQUARE_NUM,CS98_PH$REP_NUM, sep = "X")
CS07_PH$REP_ID <- paste(CS07_PH$SQUARE_NUM,CS07_PH$REP_NUM, sep = "X")
UK19_PH$REP_ID <- paste(UK19_PH$SQUARE_NUM,UK19_PH$REP_NUM, sep = "X")

PH_mod <- full_join(select(CS78_PH, REP_ID, PH1978),
                    select(CS98_PH, REP_ID, PH1998 = PHF2000)) %>%
  full_join(select(CS07_PH, REP_ID, PH2007 = PH2007_IN_WATER)) %>%
  full_join(select(UK19_PH, REP_ID, PH2019 = PH_DIW)) %>%
  pivot_longer(starts_with("PH"), values_to = "pH",
               values_drop_na = TRUE, names_prefix = "PH") %>%
  mutate(year = as.numeric(name),
         YR = as.factor(name),
         YRnm = as.integer(YR),
         SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "X"), "[",1)),
         REP_PLOT = sapply(strsplit(REP_ID, "X"), "[",2)) %>%
  left_join(select(landclass_dat, SERIES_NUM, LC07))
str(PH_mod)
summary(PH_mod)


#run the model to obtain estimates
modNat <- lme(pH ~ YR - 1, random = ~1|SERIES_NUM,
              data= PH_mod,
              correlation = corAR1(form = ~YRnm|SERIES_NUM/REP_PLOT),
              na.action = na.omit)
summary(modNat)

# check residuals 
qqp(resid(modNat), main = "pH ~ year model")
plot(modNat, main = "pH ~ year model")
plot(modNat, YR ~ resid(.))
plot(modNat, resid(., type = "p") ~ fitted(.)|YR, abline = 0,
     main = "pH ~ year model")


# checking to see if adding AVC info improves normality of residuals
# PH_mod <- PH_mod %>%
#   mutate(AggHab = recode_factor(as.factor(AVC),
#                                  "Crops/Weeds" = "Fertile",
#                                  "Fertile grassland" = "Fertile",
#                                  "Tall herb/grass" = "Fertile",
#                                  "Infertile grassland" = "Infertile",
#                                  "Lowland wooded"  = "Infertile",
#                                  "Heath/bog" = "Acidic",
#                                  "Moorland grass/mosaic" = "Acidic",
#                                  "Upland wooded" = "Acidic"))
# 
# modHab1 <- lme(pH ~ AggHab + YR - 1, random = ~1|SERIES_NUM,
#               data= PH_mod,
#               correlation = corAR1(form = ~YRnm|SERIES_NUM/REP_PLOT),
#               na.action = na.omit)
# summary(modHab1)
# 
# # check residuals 
# qqp(resid(modHab1))
# plot(modHab1)
# 
# modHab2 <- lme(pH ~ AggHab*YR - 1, random = ~1|SERIES_NUM,
#                data= PH_mod,
#                correlation = corAR1(form = ~YRnm|SERIES_NUM/REP_PLOT),
#                na.action = na.omit)
# summary(modHab2)
# 
# # check residuals 
# qqp(resid(modHab2))
# plot(modHab2)
# 
# par(mfrow=c(1,3))
# qqp(resid(modNat), main = "Year only")
# qqp(resid(modHab1), main = "AggHab + Year")
# qqp(resid(modHab2), main = "AggHab*Year")
# par(mfrow=c(1,1))

# no improvement when incorporating habitat info

## obtain yearly predictions and confidence estimates from the fitted model and store as simple table
out_dat_pred <- data.frame(Year=as.numeric(as.character(sort(unique(PH_mod$YR)))),Estimated_Value =summary(modNat)$tTable[,1],
                           Lower_est.Mod = summary(modNat)$tTable[,1]-(1.96*summary(modNat)$tTable[,2]),
                           Upper_est.Mod = summary(modNat)$tTable[,1]+(1.96*summary(modNat)$tTable[,2])
)
out_dat_pred

# pairwise comparison using emmeans
modnat.emm.s <- emmeans(modNat, "YR")
pairs(modnat.emm.s)
plot(modnat.emm.s, comparisons = TRUE)


# Bootstrapping

## obtain confidence intervals using bootstrapping (resampling within
## landclasses) as soil C is not normally distributed and not symmetrical

#to set up the bootstrapping, we create a list whereby each entry contains the
#row IDs for corresponding entried for the specific year*landclass combination
N_tab <- table(PH_mod$LC07, PH_mod$YR)
count=1
idx=list()
for(j in 1:(dim(N_tab)[1])){
  for(i in 1:(dim(N_tab)[2])){
    idx[[count]] = which(PH_mod$YR == colnames(N_tab)[i] & 
                           PH_mod$LC07 == rownames(N_tab)[j])
    count=count+1
  }
}


#create an empty matrix to store results in 
boot_ESTS=matrix(nrow=1000,ncol=4)

#we run 1000 bootstrap resamples
for(isim in 1:1000){
  
  #at each iteration we resample from the index created above - hence resampling
  #year*landclass combinations. then store the IDs
  samp_id = c()
  for(k in 1:length(idx)){
    samp_id = c(samp_id, sample(idx[[k]], length(idx[[k]]), replace=TRUE))  
  }
  
  #extract the full set of resampled row IDs and create new temporary data frame
  boot_dat = PH_mod[samp_id,]
  
  #fit model to the bootstrap sample. note that there is no longer a correlation
  #effect here. this is because we have replicates due to the resampling. this
  #doesn't matter as we are only extracting the mean values.
  boot_mod = try(lme(pH ~ YR - 1, random = ~1|SERIES_NUM,
                     data = boot_dat, na.action = na.omit))
  
  #store the yearly estimates. 
  if(class(boot_mod)!="try-error"){
    boot_ESTS[isim,] <- summary(boot_mod)$tTable[,1]
  }
  
} 

#Obtain confidence intervals from bootstrap by taking percentiles  
All_Ests <- cbind(out_dat_pred,
                  t(apply(boot_ESTS, 2, quantile, 
                          c(0.025,0.975), na.rm=TRUE)))

names(All_Ests)[5:6] <- c("Lower_Boot","Upper_Boot")


print(All_Ests)
ggplot(All_Ests, aes(x = Year)) +
  geom_line(aes(y = Estimated_Value), lty = "dotted", colour = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod), alpha = 0.1, fill = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_Boot, ymax = Upper_Boot), alpha = 0.1, fill = "#2F7ECE") +
  geom_point(aes(y = Estimated_Value), size = 3) +
  geom_linerange(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod)) +
  geom_linerange(aes(ymin = Lower_Boot, ymax = Upper_Boot), size = 2) +
  labs(y = "pH")
ggsave("pH by year model outputs bootstrap LC.png", path = "Outputs/Graphs/",
       width = 12, height = 10, units = "cm")

ggplot(All_Ests, aes(x = Year)) +
  geom_jitter(data = PH_mod, aes(x = year, y = pH), colour = "grey", alpha = 0.05, width = 1, height = 0) +
  # geom_dotplot(data = PH_mod, aes(x = year, y = pH, group = year), binaxis = "y", stackdir = "center",
  #              fill = "#2F7ECE", alpha = 0.05, width = 1, binwidth = 0.02) +
  geom_boxplot(data = PH_mod, aes(x = year, y = pH, group = year), colour = "grey20", alpha = 0.2, width = 1) +
  geom_line(aes(y = Estimated_Value), lty = "dotted", colour = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod), alpha = 0.1, fill = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_Boot, ymax = Upper_Boot), alpha = 0.1, fill = "#2F7ECE") +
  geom_point(aes(y = Estimated_Value), colour = "#2F7ECE") +
  # geom_linerange(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod)) +
  # geom_linerange(aes(ymin = Lower_Boot, ymax = Upper_Boot), size = 2) 
  labs(y = "pH")
ggsave("pH by year model bootstrap LC and data.png", path = "Outputs/Graphs/",
       width = 12, height = 10, units = "cm")



# LOI model ####
LOI <- full_join(select(CS78_PH, REP_ID, LOI1978),
                 select(CS98_PH, REP_ID, LOI1998 = LOI2000)) %>%
  full_join(select(CS07_PH, REP_ID, LOI2007)) %>%
  full_join(select(UK19_PH, REP_ID, LOI2019 = LOI)) %>%
  pivot_longer(starts_with("LOI"), values_to = "LOI",
               values_drop_na = TRUE, names_prefix = "LOI") %>%
  mutate(year = as.numeric(name),
         YR = as.factor(name),
         YRnm = as.integer(YR),
         SERIES_NUM = sapply(strsplit(REP_ID, "X"), "[",1),
         REP_PLOT = sapply(strsplit(REP_ID, "X"), "[",2),
         LOI_pr = 0.01*LOI,
         CARBO_LOI = 5.5*LOI)
str(LOI)
summary(LOI)

#run the model to obtain estimates
modLOI <- lme(LOI ~ YR - 1, random = ~1|SERIES_NUM,
              data = LOI,
              correlation = corAR1(form = ~YRnm|SERIES_NUM/REP_PLOT),
              na.action = na.omit)
qqp(resid(modLOI), main = "LOI ~ year model")
plot(modLOI, main = "LOI ~ year model")



# pairwise comparison using emmeans
modloi.emm.s <- emmeans(modLOI, "YR")
pairs(modloi.emm.s)
plot(modloi.emm.s, comparisons = TRUE)
# modLOI <- lme(log(LOI) ~ YR - 1, random = ~1|SERIES_NUM,
#               data = LOI,
#               correlation = corAR1(form = ~YRnm|SERIES_NUM/REP_PLOT),
#               na.action = na.omit)
# qqp(resid(modLOI), main = "log(LOI) ~ year model")
# plot(modLOI, main = "log(LOI) ~ year model")
# 
library(glmmTMB)
modLOI <- glmmTMB(LOI_pr ~ YR - 1 + (1|SERIES_NUM) +
                    ar1(YRnm + 0|REP_ID),
                  data= LOI, family = beta_family(),
                  na.action = na.omit)
summary(modLOI)
# 
# 
# modLOI_noar <- glmmTMB(LOI_pr ~ YR - 1 + (1|SERIES_NUM),
#                   data= LOI, family = beta_family(),
#                   na.action = na.omit)
# 
# # check residuals using DHARMa 
resLOI_noar <- simulateResiduals(modLOI_noar)
plot(resLOI_noar)

# Bootstrapping

## obtain confidence intervals using bootstrapping (resampling within
## landclasses) as soil C is not normally distributed and not symmetrical

#to set up the bootstrapping, we create a list whereby each entry contains the
#row IDs for corresponding entried for the specific year*landclass combination
LOI <- LOI %>% mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "X"), "[", 1))) %>%
  left_join(select(landclass_dat, SERIES_NUM, LC07))
summary(LOI)

N_tab <- table(LOI$LC07, LOI$YR)
count=1
idx=list()
for(j in 1:(dim(N_tab)[1])){
  for(i in 1:(dim(N_tab)[2])){
    idx[[count]] = which(LOI$YR == colnames(N_tab)[i] & 
                           LOI$LC07 == rownames(N_tab)[j])
    count=count+1
  }
}


#create an empty matrix to store results in 
boot_ESTS=matrix(nrow=1000,ncol=4)

#we run 1000 bootstrap resamples
for(isim in 1:1000){
  
  #at each iteration we resample from the index created above - hence resampling
  #year*landclass combinations. then store the IDs
  samp_id = c()
  for(k in 1:length(idx)){
    samp_id = c(samp_id, sample(idx[[k]], length(idx[[k]]), replace=TRUE))  
  }
  
  #extract the full set of resampled row IDs and create new temporary data frame
  boot_dat = LOI[samp_id,]
  
  #fit model to the bootstrap sample. note that there is no longer a correlation
  #effect here. this is because we have replicates due to the resampling. this
  #doesn't matter as we are only extracting the mean values.
  boot_mod = try(lme(LOI ~ YR - 1, random = ~1|SERIES_NUM,
                     data = boot_dat, na.action = na.omit))
  
  #store the yearly estimates. 
  if(class(boot_mod)!="try-error"){
    boot_ESTS[isim,] <- summary(boot_mod)$tTable[,1]
  }
  
} 

out_dat_pred <- data.frame(Year=as.numeric(as.character(sort(unique(LOI$YR)))),Estimated_Value =summary(modLOI)$tTable[,1],
                           Lower_est.Mod = summary(modLOI)$tTable[,1]-(1.96*summary(modLOI)$tTable[,2]),
                           Upper_est.Mod = summary(modLOI)$tTable[,1]+(1.96*summary(modLOI)$tTable[,2])
)
out_dat_pred

#Obtain confidence intervals from bootstrap by taking percentiles  
All_Ests <- cbind(out_dat_pred,
                  t(apply(boot_ESTS, 2, quantile, 
                          c(0.025,0.975), na.rm=TRUE)))

names(All_Ests)[5:6] <- c("Lower_Boot","Upper_Boot")


print(All_Ests)
ggplot(All_Ests, aes(x = Year)) +
  geom_line(aes(y = Estimated_Value), lty = "dotted", colour = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod), alpha = 0.1, fill = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_Boot, ymax = Upper_Boot), alpha = 0.1, fill = "#2F7ECE") +
  geom_point(aes(y = Estimated_Value), size = 3) +
  geom_linerange(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod)) +
  geom_linerange(aes(ymin = Lower_Boot, ymax = Upper_Boot), size = 2) +
  labs(y = "LOI (%)")
ggsave("LOI by year model outputs bootstrap LC.png", path = "Outputs/Graphs/",
       width = 12, height = 10, units = "cm")

ggplot(All_Ests, aes(x = Year)) +
  geom_jitter(data = LOI, aes(x = year, y = LOI), colour = "grey", alpha = 0.05, width = 1, height = 0) +
  geom_boxplot(data = LOI, aes(x = year, y = LOI, group = year), colour = "grey20", alpha = 0.2, width = 1) +
  geom_line(aes(y = Estimated_Value), lty = "dotted", colour = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod), alpha = 0.1, fill = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_Boot, ymax = Upper_Boot), alpha = 0.1, fill = "#2F7ECE") +
  geom_point(aes(y = Estimated_Value), colour = "#2F7ECE") +
  labs(y = "LOI (%)")
ggsave("LOI by year model bootstrap LC and data.png", path = "Outputs/Graphs/",
       width = 12, height = 10, units = "cm")


# LOI log ####
modLOI <- lme(log(LOI) ~ YR - 1, random = ~1|SERIES_NUM,
              data = LOI,
              correlation = corAR1(form = ~YRnm|SERIES_NUM/REP_PLOT),
              na.action = na.omit)
qqp(resid(modLOI), main = "log(LOI) ~ year model")
plot(modLOI, main = "log(LOI) ~ year model")

# pairwise comparison using emmeans
modloi.log.emm.s <- emmeans(modLOI, "YR")
pairs(modloi.log.emm.s)
plot(modloi.log.emm.s, comparisons = TRUE)

# Bootstrapping

## obtain confidence intervals using bootstrapping (resampling within
## landclasses) as soil C is not normally distributed and not symmetrical

#to set up the bootstrapping, we create a list whereby each entry contains the
#row IDs for corresponding entried for the specific year*landclass combination

N_tab <- table(LOI$LC07, LOI$YR)
count=1
idx=list()
for(j in 1:(dim(N_tab)[1])){
  for(i in 1:(dim(N_tab)[2])){
    idx[[count]] = which(LOI$YR == colnames(N_tab)[i] & 
                           LOI$LC07 == rownames(N_tab)[j])
    count=count+1
  }
}


#create an empty matrix to store results in 
boot_ESTS=matrix(nrow=1000,ncol=4)

#we run 1000 bootstrap resamples
for(isim in 1:1000){
  
  #at each iteration we resample from the index created above - hence resampling
  #year*landclass combinations. then store the IDs
  samp_id = c()
  for(k in 1:length(idx)){
    samp_id = c(samp_id, sample(idx[[k]], length(idx[[k]]), replace=TRUE))  
  }
  #extract the full set of resampled row IDs and create new temporary data frame
  boot_dat = LOI[samp_id,]
  
  #fit model to the bootstrap sample. note that there is no longer a correlation
  #effect here. this is because we have replicates due to the resampling. this
  #doesn't matter as we are only extracting the mean values.
  boot_mod = try(lme(log(LOI) ~ YR - 1, random = ~1|SERIES_NUM,
                     data = boot_dat, na.action = na.omit))
  
  #store the yearly estimates. 
  if(class(boot_mod)!="try-error"){
    boot_ESTS[isim,] <- summary(boot_mod)$tTable[,1]
  }
  
} 

out_dat_pred_log <- data.frame(Year=as.numeric(as.character(sort(unique(LOI$YR)))),Estimated_Value =summary(modLOI)$tTable[,1],
                               Lower_est.Mod = summary(modLOI)$tTable[,1]-(1.96*summary(modLOI)$tTable[,2]),
                               Upper_est.Mod = summary(modLOI)$tTable[,1]+(1.96*summary(modLOI)$tTable[,2])
)
out_dat_pred_log

#Obtain confidence intervals from bootstrap by taking percentiles  
All_Ests_log <- cbind(out_dat_pred_log,
                      t(apply(boot_ESTS, 2, quantile, 
                              c(0.025,0.975), na.rm=TRUE)))

names(All_Ests_log)[5:6] <- c("Lower_Boot","Upper_Boot")

All_Ests_log_tr <- as.data.frame(
  cbind(Year = All_Ests_log[,"Year"],
        apply(All_Ests_log[,2:6],1:2, exp)))

print(All_Ests_log_tr)
ggplot(All_Ests_log_tr, aes(x = Year)) +
  geom_line(aes(y = Estimated_Value), lty = "dotted", colour = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod), alpha = 0.1, fill = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_Boot, ymax = Upper_Boot), alpha = 0.1, fill = "#2F7ECE") +
  geom_point(aes(y = Estimated_Value), size = 3) +
  geom_linerange(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod)) +
  geom_linerange(aes(ymin = Lower_Boot, ymax = Upper_Boot), size = 2) +
  labs(y = "LOI (%)")
ggsave("LOI by year log model outputs bootstrap LC.png", path = "Outputs/Graphs/",
       width = 12, height = 10, units = "cm")

LOI_mean <- LOI %>%
  group_by(year) %>%
  summarise(mean_LOI = mean(LOI),
            se_LOI = sd(LOI)/sqrt(length(LOI))) %>%
  mutate(upper_mean = mean_LOI + se_LOI,
         lower_mean = mean_LOI - se_LOI)

ggplot(All_Ests_log_tr, aes(x = Year)) +
  geom_jitter(data = LOI, aes(x = year, y = LOI), colour = "grey", alpha = 0.05, width = 1, height = 0) +
  geom_boxplot(data = LOI, aes(x = year, y = LOI, group = year), colour = "grey20", alpha = 0.2, width = 1) +
  # geom_pointrange(data = LOI_mean, aes(x = year, y = mean_LOI,
  #                                 ymax = upper_mean, ymin = lower_mean), size = 0.1) +
  
  # log
  geom_line(aes(y = Estimated_Value), lty = "dotted", colour = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod), alpha = 0.1, fill = "#2F7ECE") +
  geom_ribbon(aes(ymin = Lower_Boot, ymax = Upper_Boot), alpha = 0.1, fill = "#2F7ECE") +
  geom_point(aes(y = Estimated_Value), colour = "#2F7ECE") +
  
  # non log
  geom_line(data = All_Ests, aes(y = Estimated_Value), lty = "dotted", colour = "#4DA43A") +
  geom_ribbon(data = All_Ests, aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod), alpha = 0.1, fill = "#4DA43A") +
  geom_ribbon(data = All_Ests, aes(ymin = Lower_Boot, ymax = Upper_Boot), alpha = 0.1, fill = "#4DA43A") +
  geom_point(data = All_Ests, aes(y = Estimated_Value), colour = "#4DA43A") +
  
  labs(y = "LOI (%)")
ggsave("LOI by year log blue nolog green model bootstrap LC and data.png", path = "Outputs/Graphs/",
       width = 12, height = 10, units = "cm")

# Bayesian version ####
library(brms)
get_prior(pH ~ YR -1 + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID), data = PH_mod)

mod_pr <- prior(normal(5.5,1), class = "b")

brm_pH <- brm(pH ~ YR -1 + (1|SERIES_NUM) + 
                ar(time = YRnm, gr = REP_ID), data = PH_mod,
              prior = mod_pr, cores = 4, iter = 4000)


get_prior(LOI_pr ~ YR - 1 + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID, cov = TRUE), 
          data = LOI, family = "beta")

mod_pr <- c(prior(normal(0.3,0.1), class = "b"),
            prior(normal(0.2,0.1), class = "ar"),
            prior(gamma(0.1,0.1), class = "phi"),
            prior(student_t(3, 0, 1), class = "sd"))
brm_LOI <- brm(LOI_pr ~ YR -1 + (1|SERIES_NUM) + 
                 ar(time = YRnm, gr = REP_ID, cov = TRUE), 
               data = LOI, family = "beta",
               prior = mod_pr, cores = 2, chains = 2, iter = 2000,
               sample_prior = "only")
plot(brm_LOI)
pp_check(brm_LOI)




# AVC
LOI <- LOI %>%
  mutate(Acidity = recode_factor(as.factor(AVC),
                                 "Crops/Weeds" = "Neutral",
                                 "Fertile grassland" = "Neutral",
                                 "Tall herb/grass" = "Neutral",
                                 "Infertile grassland" = "Neutral",
                                 "Lowland wooded"  = "Neutral",
                                 "Heath/bog" = "Acidic",
                                 "Moorland grass/mosaic" = "Neutral",
                                 "Upland wooded" = "Neutral")) %>%
  mutate(Acidity =tidyr::replace_na(Acidity, "Neutral"))

ggplot(LOI, aes(x = LOI)) + geom_histogram() + facet_wrap(~Acidity)
get_prior(LOI_pr ~ YR - 1 + Acidity + (1|SERIES_NUM) + 
            ar(time = YRnm, gr = REP_ID, cov = TRUE), 
          data = LOI, family = "beta")

mod_pr <- c(prior(normal(-2,0.5), class = "b", coef = "YR1978"),
            prior(normal(-2,0.5), class = "b", coef = "YR1998"),
            prior(normal(-2,0.5), class = "b", coef = "YR2007"),
            prior(normal(-2,0.5), class = "b", coef = "YR2019"),
            prior(normal(3,0.5), class = "b", coef = "AcidityAcidic"),
            prior(normal(0.2,0.2), class = "ar"),
            prior(gamma(100,0.2), class = "phi"),
            prior(student_t(3, 0, 1), class = "sd"),
            prior(student_t(3, 0, 0.5), class = "sderr"))
brm_LOI2 <- brm(LOI_pr ~ YR -1 + Acidity + (1|SERIES_NUM) + 
                  ar(time = YRnm, gr = REP_ID, cov = TRUE), 
                data = LOI, family = "beta",
                prior = mod_pr, cores = 2, chains = 4, iter = 4000)
saveRDS(brm_LOI2, "Outputs/Models/brm_LOI2_070720_bfmiwarning.rds")

# checking energy correlations
test <- get_sampler_params(brm_LOI2$fit, inc_warmup = FALSE)
test2 <- do.call(rbind, test)
brm_samps <- as.matrix(brm_LOI2)

for(i in 1:ncol(brm_samps)){
  cor2 <- cor(brm_samps[,i],test2[,"energy__"])
  if(abs(cor2)>=0.1)
    print(paste(colnames(brm_samps)[i],cor2))
}
# [1] "phi -0.891683685248258"
# [1] "lp__ -0.963763534407218"
# [1] "sderr -0.222745655770828"

# lp is supposed to be correlated, the issue is phi (maybe sderr?)
test3 <- cbind(brm_samps[,c("phi","sderr")],test2[,"energy__"])

pairs(test3)

# with data
mod_pr <- c(prior(normal(-1,0.3), class = "b"),
            prior(normal(0.2,0.2), class = "ar"),
            prior(gamma(8,12), class = "phi"),
            prior(student_t(3, 0, 1), class = "sd"),
            prior(student_t(3, 0, 0.5), class = "sderr"))
brm_LOI <- brm(LOI_pr ~ YR -1 + (1|SERIES_NUM) + 
                 ar(time = YRnm, gr = REP_ID, cov = TRUE), 
               data = LOI, family = "beta",
               prior = mod_pr, chains = 4, cores = 4, iter = 2000)

# model checks
summary(brm_LOI)
plot(brm_LOI, ask = FALSE)
pp_check(brm_LOI)
pp_check(brm_LOI, type = "stat")

plot(conditional_effects(brm_LOI))

modLOI <- lme(LOI ~ YR - 1, random = ~1|SERIES_NUM,
              data = LOI,
              correlation = corAR1(form = ~YRnm|SERIES_NUM/REP_PLOT),
              na.action = na.omit)

out_dat_pred <- data.frame(Year=as.numeric(as.character(sort(unique(LOI$YR)))),Estimated_Value =summary(modLOI)$tTable[,1],
                           Lower_est.Mod = summary(modLOI)$tTable[,1]-(1.96*summary(modLOI)$tTable[,2]),
                           Upper_est.Mod = summary(modLOI)$tTable[,1]+(1.96*summary(modLOI)$tTable[,2])
)
out_dat_pred

# bayes estimat
brm_fixef <- fixef(brm_LOI, robust = FALSE)
brm_fixef <- apply(brm_fixef, 1:2, function(x) 100*exp(x)/(1 + exp(x)))

brm_fixef_med <- fixef(brm_LOI, robust = TRUE)
brm_fixef_med <- apply(brm_fixef_med, 1:2, function(x) 100*exp(x)/(1 + exp(x)))
colnames(brm_fixef_med) <- paste0(colnames(brm_fixef_med), "_med") 

colnames(out_dat_pred_log) <- paste0(colnames(out_dat_pred_log),"_log")

out_dat_pred_log <- exp(out_dat_pred_log)
#Obtain confidence intervals from bootstrap by taking percentiles  
All_Ests <- cbind(out_dat_pred,
                  out_dat_pred_log[,-1],
                  brm_fixef[,-2],
                  brm_fixef_med[,-2])

All_Ests

ggplot(All_Ests, aes(x = Year)) +
  geom_jitter(data = LOI, aes(x = year, y = LOI), colour = "grey", alpha = 0.05, width = 1, height = 0) +
  geom_boxplot(data = LOI, aes(x = year, y = LOI, group = year), 
               colour = "grey20", alpha = 0.2, width = 1,
               outlier.shape = NA) +
  geom_pointrange(data = LOI_mean, aes(x = year, y = mean_LOI,
                                       ymax = upper_mean, ymin = lower_mean), size = 0.1) +
  
  # Bayes - blue
  geom_line(aes(y = Estimate), lty = "dotted", colour = "#2F7ECE") +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.1, fill = "#2F7ECE") +
  geom_point(aes(y = Estimate), colour = "#2F7ECE") +
  
  # # Bayes - purple
  # geom_line(aes(y = Estimate_med), lty = "dotted", colour = "purple") +
  # geom_ribbon(aes(ymin = Q2.5_med, ymax = Q97.5_med), alpha = 0.1, fill = "purple") +
  # geom_point(aes(y = Estimate_med), colour = "purple") +
  # 
  # nlme - green
  geom_line(aes(y = Estimated_Value), lty = "dotted", colour = "#4DA43A") +
  geom_ribbon(aes(ymin = Lower_est.Mod, ymax = Upper_est.Mod), alpha = 0.1, fill = "#4DA43A") +
  geom_point(aes(y = Estimated_Value), colour = "#4DA43A") +
  
  # nlme log - pink
  geom_line(aes(y = Estimated_Value_log), lty = "dotted", colour = "pink") +
  geom_ribbon(aes(ymin = Lower_est.Mod_log, ymax = Upper_est.Mod_log), 
              alpha = 0.1, fill = "pink") +
  geom_point(aes(y = Estimated_Value_log), colour = "pink") +
  
  labs(y = "LOI (%)")

ggplot(LOI, aes(x =LOI)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~year, nrow = 4, scales = "free_y") + 
  scale_x_continuous(expand=c(0,0)) + 
  labs(x = "LOI %")
ggsave("LOI histogram by year.png", path = "Outputs/Graphs/",
       width = 10, height = 15, units = "cm")


# Mixture model
mix <- mixture(Beta, Beta)
mod_pr <- c(
  prior(normal(-2,1), class = "Intercept", dpar = "mu1"),
  prior(normal(2,1), class = "Intercept", dpar = "mu2"),
  prior(normal(0,.5), class = "b", dpar = "mu1"),
  prior(normal(0,.5), class = "b", dpar = "mu2"),
  prior(gamma(0.01,0.01), class = "phi1"),
  prior(gamma(0.01,0.01), class = "phi2")
)

brm_LOI_mix <- brm(bf(LOI_pr ~ YR), data = LOI, family = mix,
                   prior = mod_pr, chains = 4, cores = 4)
summary(brm_LOI_mix)

plot(brm_LOI_mix)

pp_check(brm_LOI_mix, nsamples = 30) 
# not the best recovery of the data - possibly a third peak at ~ 30% is messing it up

pp_check(brm_LOI_mix, type = "stat_2d")

plot(conditional_effects(brm_LOI_mix))

# mixture model - spatiotemporal
mod_pr <- c(
  prior(normal(-2,1), class = "Intercept", dpar = "mu1"),
  prior(normal(2,1), class = "Intercept", dpar = "mu2"),
  prior(normal(0,.5), class = "b", dpar = "mu1"),
  prior(normal(0,.5), class = "b", dpar = "mu2"),
  prior(gamma(0.01,0.01), class = "phi1"),
  prior(gamma(0.01,0.01), class = "phi2"),
  # prior(normal(0.2,0.2), class = "ar"),
  prior(student_t(3, 0, 0.2), class = "sd", dpar = "mu1"),
  prior(student_t(3, 0, 0.2), class = "sd", dpar = "mu2")
  # prior(student_t(3, 0, 0.5), class = "sderr")
)

brm_LOI_mixc <- brm(bf(LOI_pr ~ YR + (1|SERIES_NUM)), 
                    data = LOI, family = mix,
                    prior = mod_pr, chains = 4, cores = 4)

# mixture model - 3 groups
mix3 <- mixture(Beta, Beta, Beta)
mod_pr3 <- c(
  prior(normal(-2.5,0.5), class = "Intercept", dpar = "mu1"),
  prior(normal(-0.5,0.5), class = "Intercept", dpar = "mu2"),
  prior(normal(2.5,0.5), class = "Intercept", dpar = "mu3"),
  prior(normal(0,.5), class = "b", dpar = "mu1"),
  prior(normal(0,.5), class = "b", dpar = "mu2"),
  prior(normal(0,.5), class = "b", dpar = "mu3"),
  prior(gamma(0.01,0.01), class = "phi1"),
  prior(gamma(0.01,0.01), class = "phi2"),
  prior(gamma(0.01,0.01), class = "phi3")
)

brm_LOI_mix3 <- brm(bf(LOI_pr ~ YR), data = LOI, family = mix3,
                    prior = mod_pr3, chains = 4, cores = 4)
summary(brm_LOI_mix3)

# mixture model - 3 groups + spatial
mix3 <- mixture(Beta, Beta, Beta)
mod_pr3 <- c(
  prior(normal(-2.5,0.5), class = "Intercept", dpar = "mu1"),
  prior(normal(-0.5,0.5), class = "Intercept", dpar = "mu2"),
  prior(normal(2.5,0.5), class = "Intercept", dpar = "mu3"),
  prior(normal(0,.5), class = "b", dpar = "mu1"),
  prior(normal(0,.5), class = "b", dpar = "mu2"),
  prior(normal(0,.5), class = "b", dpar = "mu3"),
  prior(gamma(0.01,0.01), class = "phi1"),
  prior(gamma(0.01,0.01), class = "phi2"),
  prior(gamma(0.01,0.01), class = "phi3"),
  prior(student_t(3, 0, 0.2), class = "sd", dpar = "mu1"),
  prior(student_t(3, 0, 0.2), class = "sd", dpar = "mu2"),
  prior(student_t(3, 0, 0.2), class = "sd", dpar = "mu3")
)

brm_LOI_mix3b <- brm(bf(LOI_pr ~ YR + (1|SERIES_NUM)), 
                     data = LOI, family = mix3,
                    prior = mod_pr3, chains = 4, cores = 4)
summary(brm_LOI_mix3b)

# manually edited code to use only one sd for all peaks
brm_data <- standata(brm_LOI_mix3b)

library(rstan)
LOI_mix3bv2 <- stan(
  "beta_3mix_sq.stan", data = brm_data, chains = 4, cores = 4
)
summary(LOI_mix3bv2)