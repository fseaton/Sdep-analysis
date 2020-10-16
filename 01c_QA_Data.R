# Script for QA analysis for Ellenberg R and soil pH
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())

# Ellenberg QA data ####
# Only doing this for Ellenberg R right now
# for small plot calculations can also include Y and U plots (both 4m2)
str(CSVEG_QA)

CS16_VEG_QA2 <- CS16_VEG_QA %>%
  mutate(REP_ID = substring(Quadrat,2)) %>%
  mutate(PLOT_TYPE = gsub("[0-9]","",REP_ID)) %>%
  filter(PLOT_TYPE %in% c("X","Y","U")) %>%
  left_join(select(VEGETATION_PLOTS_20161819, REP_ID = REP_ID_NEW,
                   CS_REP_ID = REP_ID)) %>%
  select(-REP_ID) %>%
  select(REP_ID = CS_REP_ID, PLOT_TYPE, BRC_NUMBER = Names, COVER = Cover,
         Surveyor = SURVEY) %>%
  mutate(Year = 2019,
         Surveyor = recode(Surveyor, 
                           "SV" = "CS"))

CS_VEG_QA2 <- CS_VEG_QA %>% 
  mutate(REP_ID = paste0(Square, Plot_Type, Plot_number)) %>%
  filter(Plot_Type %in% c("X","Y","U")) %>%
  select(REP_ID, Year, PLOT_TYPE = Plot_Type,
         BRC_NUMBER, COVER = Value, Surveyor)
  

X_Ell_QA_whole <- CS19_VEG_QA %>%
  filter(Surveyor == "QA") %>%
  select(REP_ID, PLOT_TYPE, BRC_NUMBER, COVER = TOTAL_NUM) %>%
  mutate(Year = 2019, Surveyor = "QA") %>%
  full_join(CS16_VEG_QA2) %>%
  full_join(CS_VEG_QA2) %>%
  filter(PLOT_TYPE == "X") %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   EBERGR)) %>%
  mutate(across(EBERGR, na_if, y = 0)) %>% # set 0 values to NA
  filter(!is.na(EBERGR)) %>%
  group_by(REP_ID, Year, Surveyor) %>%
  summarise(across(EBERGR, 
                   .fns = list(R_UW = ~mean(.x, na.rm = TRUE),
                               R_W = ~weighted.mean(.x, w = COVER, na.rm = TRUE)))) %>%
  rename_with(~gsub("EBERGR","WH",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA))

X_Ell_QA_inner <- CS19_VEG_QA %>%
  filter(Surveyor == "QA" & NEST_LEVEL < 2) %>%
  select(REP_ID, PLOT_TYPE, BRC_NUMBER, COVER = FIRST_COVER) %>%
  mutate(Year = 2019, Surveyor = "QA") %>%
  full_join(filter(CS16_VEG_QA2, PLOT_TYPE %in% c("Y","U"))) %>%
  full_join(filter(CS_VEG_QA2, PLOT_TYPE %in% c("Y","U"))) %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   EBERGR)) %>%
  mutate(across(EBERGR, na_if, y = 0)) %>% # set 0 values to NA
  filter(!is.na(EBERGR)) %>%
  group_by(REP_ID, Year, Surveyor) %>%
  summarise(across(EBERGR, 
                   .fns = list(R_UW = ~mean(.x, na.rm = TRUE),
                               R_W = ~weighted.mean(.x, w = COVER, na.rm = TRUE)))) %>%
  rename_with(~gsub("EBERGR","SM",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA))

X_Ell_QA <- full_join(X_Ell_QA_inner, X_Ell_QA_whole) %>%
  left_join(CS_REP_ID_LONG) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID)) %>%
  select(-REPEAT_PLOT_ID) %>%
  left_join(select(X_Ell, REP_ID, Year,
                   SM_R_W, SM_R_UW,WH_R_W,WH_R_UW) %>%
              filter(Year == 2019) %>%
              mutate(Surveyor = "CS")) %>%
  pivot_longer(contains("_R_"), 
               names_to = c("PlotSize","Ellenberg","Weighting"),
               names_sep = "_", values_to = "score") %>%
  filter(!is.na(score)) %>%
  pivot_wider(names_from = "Surveyor", values_from = "score") %>%
  mutate(Diff = QA - CS)


ggplot(X_Ell_QA, aes(x = CS, y = QA, colour = Year)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(PlotSize ~ Weighting)

ggplot(X_Ell_QA, aes(x = Diff)) + 
  geom_histogram() +
  facet_wrap(~Year + PlotSize + Weighting)


janitor::tabyl(X_Ell_QA, Year, PlotSize, Weighting)
# $UW
# Year SM WH
# 1990 39 44
# 1998 51 39
# 2007 63 52
# 2019 43 43
# 
# $W
# Year SM WH
# 1990 39 44
# 1998 51 39
# 2007 63 52
# 2019 43 43

X_Ell_QA %>% 
  # mutate(Year = ifelse(Year == 2016, 2019, Year)) %>%
  group_by(Year, PlotSize, Weighting) %>%
  na.omit() %>%
  summarise(mu = MASS::fitdistr(Diff,"t", 
                                start = list(m = 0, s = 0.2, df = 5),
                                lower=c(-1, 0.001,1))$estimate[1],
            sd = MASS::fitdistr(Diff,"t", 
                                start = list(m = 0, s = 0.2, df = 5),
                                lower=c(-1, 0.001,1))$estimate[2],
            df = MASS::fitdistr(Diff,"t", 
                                start = list(m = 0, s = 0.2, df = 5),
                                lower=c(-1, 0.001,1))$estimate[3])
# PlotSize Weighting      mu    sd    df
# <chr>    <chr>       <dbl> <dbl> <dbl>
# 1 SM       UW         0.0452 0.256  4.07
# 2 SM       W          0.0373 0.378  6.64
# 3 WH       UW        -0.0388 0.128  5.50
# 4 WH       W         -0.0714 0.298 29.9 

X_Ell_QA_hab <- left_join(X_Ell_QA, BH_comb_nodupes)
table(X_Ell_QA_hab$BH_DESC)



X_Ell_QA_hab <- filter(X_Ell_QA_hab,
                       BH %in% c(1,5,6,7,8,9,10,11,12)) %>%
  mutate(Improved = ifelse(BH %in% c(5,6),1,0))

janitor::tabyl(X_Ell_QA_hab, Year, PlotSize, Improved)
X_Ell_QA_hab %>% 
  group_by(Year, PlotSize, Weighting) %>%
  na.omit() %>%
  summarise(mu = MASS::fitdistr(Diff,"t", 
                                start = list(m = 0, s = 0.2, df = 5),
                                lower=c(-1, 0.001,1))$estimate[1],
            sd = MASS::fitdistr(Diff,"t", 
                                start = list(m = 0, s = 0.2, df = 5),
                                lower=c(-1, 0.001,1))$estimate[2],
            df = MASS::fitdistr(Diff,"t", 
                                start = list(m = 0, s = 0.2, df = 5),
                                lower=c(-1, 0.001,1))$estimate[3])



Ell_9098 <- brms::rstudent_t(100000,6.95,0.103,0.169) -
  brms::rstudent_t(100000,2.58,0.00775,0.121)
MASS::fitdistr(Ell_9098, "t")
# m              s              df     
# 0.0955785286   0.2241806089   4.4623073358 
# (0.0008300429) (0.0008397081) (0.0605551027)
Ell_9807 <- brms::rstudent_t(100000,4.08,0.0111,0.164) - 
  brms::rstudent_t(100000,6.95,0.103,0.169)
MASS::fitdistr(Ell_9807, "t")
# m               s               df      
# -0.0914398076    0.2526516756    6.5811814098 
# ( 0.0008985194) ( 0.0009267179) ( 0.1258995813)
Ell_0719 <- brms::rstudent_t(100000,5.25,0.0169,0.150) - 
  brms::rstudent_t(100000,4.08,0.0111,0.164) 
MASS::fitdistr(Ell_0719, "t")
# m              s              df     
# 0.0051229646   0.2419197912   5.9533278679 
# (0.0008683201) (0.0009036674) (0.1059707380)

# check normal distribution
X_Ell_whole_QA %>% 
  mutate(Year = ifelse(Year == 2016, 2019, Year)) %>%
  group_by(Year) %>%
  na.omit() %>%
  summarise(mu = MASS::fitdistr(WH_UW_diff,"normal")$estimate[1],
            sd = MASS::fitdistr(WH_UW_diff,"normal")$estimate[2])
# Year       mu    sd
# <dbl>    <dbl> <dbl>
# 1  1990 -0.00476 0.203
# 2  1998  0.0936  0.198
# 3  2007  0.0248  0.236
# 4  2019  0.0115  0.185

ELL_QA_WH_NORM <- X_Ell_whole_QA %>% 
  mutate(Year = ifelse(Year == 2016, 2019, Year)) %>% 
  rbind(mutate(X_Ell_whole_QA, Year = 1978)) %>%
  group_by(Year) %>%
  na.omit() %>%
  summarise(sd = MASS::fitdistr(WH_UW_diff,"normal")$estimate[2])


# Check by habitat
X_Ell_whole_QA_BH <- do.call(rbind, list(
  X_Ell_whole_QA %>% filter(Year == 1990 & !is.na(WH_UW_diff)) %>%
    select(REP_ID, Year, WH_UW_diff) %>%
    left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y90)),
  X_Ell_whole_QA %>% filter(Year == 1998 & !is.na(WH_UW_diff)) %>%
    select(REP_ID, Year, WH_UW_diff) %>%
    left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y9899)),
  X_Ell_whole_QA %>% filter(Year == 2007 & !is.na(WH_UW_diff)) %>%
    select(REP_ID, Year, WH_UW_diff) %>%
    left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y07)),
  X_Ell_whole_QA %>% filter(Year %in% c(2016,2019) & !is.na(WH_UW_diff)) %>%
    select(REP_ID, Year, WH_UW_diff) %>%
    left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y19))
)) %>%
  select(-REP_ID) %>%
  left_join(BH_comb, by = c("REPEAT_PLOT_ID" = "REP_ID","Year")) %>%
  mutate(Habitat = ifelse(BH %in% c(1,2), 
                          "Woodland",
                          ifelse(BH %in% c(4,5,6),
                                 "Improved", 
                                 ifelse(BH %in% c(7,8,9,10,11,12,15,16), 
                                        "Habitat","Other"))))



ggplot(X_Ell_whole_QA_BH, aes(x = BH_DESC, y = WH_UW_diff)) +
  geom_jitter(height = 0, colour = "grey") +
  geom_boxplot(fill= NA, outlier.shape = NA)
ggplot(X_Ell_whole_QA_BH, aes(x = Habitat, y = WH_UW_diff)) +
  geom_jitter(height = 0, colour = "grey") +
  geom_boxplot(fill= NA, outlier.shape = NA)


X_Ell_whole_QA_BH %>% 
  mutate(Year = ifelse(Year == 2016, 2019, Year)) %>%
  group_by(Year, Habitat) %>%
  na.omit() %>%
  summarise(mu = MASS::fitdistr(WH_UW_diff,"normal")$estimate[1],
            sd = MASS::fitdistr(WH_UW_diff,"normal")$estimate[2])
# Year Habitat        mu      sd
# <dbl> <chr>       <dbl>   <dbl>
#   1  1990 Habitat   0.0264   0.0997
# 2  1990 Improved -0.0173   0.247 
# 3  1990 Other     0.0303  NA     
# 4  1990 Woodland  0.00257  0.0442
# 5  1998 Habitat   0.172    0.128 
# 6  1998 Improved  0.0379   0.217 
# 7  1998 Other     0.317   NA     
# 8  1998 Woodland  0.195    0.0358
# 9  2007 Habitat  -0.0148   0.177 
# 10  2007 Improved  0.0179   0.190 
# 11  2007 Other     0.0891  NA     
# 12  2007 Woodland  0.338    0.519 
# 13  2019 Habitat   0.118    0.155 
# 14  2019 Improved -0.0143   0.190 
# 15  2019 Woodland  0.0444  NA     
X_Ell_whole_QA_BH %>% 
  mutate(Year = ifelse(Year == 2016, 2019, Year)) %>%
  filter(Habitat != "Woodland") %>%
  filter(Year != 1990) %>%
  group_by(Year, Habitat) %>%
  na.omit() %>%
  summarise(mu = MASS::fitdistr(WH_UW_diff,"t", start = list(m = 0, s = 0.5, df = 3))$estimate[1],
            sd = MASS::fitdistr(WH_UW_diff,"t", start = list(m = 0, s = 0.5, df = 3))$estimate[2],
            df = MASS::fitdistr(WH_UW_diff,"t", start = list(m = 0, s = 0.5, df = 3))$estimate[3])




# pH QA data ####
# 1978 to 1998
str(CS78_PH_QA)
CS78_PH_QA <- mutate(CS78_PH_QA, 
                     across(.fns = ~round(.x/0.005)*0.005))

ggplot(CS78_PH_QA, aes(x = PH_1971, y = PH_2000)) +
  geom_point() + 
  geom_smooth(method = "lm")
summary(lm(PH_2000 ~ PH_1971, CS78_PH_QA))

ph78_diff <- CS78_PH_QA$PH_2000 - CS78_PH_QA$PH_1971 
MASS::fitdistr(ph78_diff, "t")
# m             s            df     
# -0.29246641    0.35713021    4.65432249 
# ( 0.03176037) ( 0.03581546) ( 1.81106541)
hist(ph78_diff)

# compare visually fits of normal and student T distribution
h <- hist(ph78_diff, breaks = 40)
xfit<-seq(min(ph78_diff),max(ph78_diff),length=40)
yfit<-dnorm(xfit,mean=mean(ph78_diff),sd=sd(ph78_diff))
yfit <- yfit*diff(h$mids[1:2])*length(ph78_diff)
lines(xfit, yfit, col="blue", lwd=2)
yfit<-brms::dstudent_t(xfit,mu = -0.29246641, sigma = 0.35713021, df = 4.65432249)
yfit <- yfit*diff(h$mids[1:2])*length(ph78_diff)
lines(xfit, yfit, col="red", lwd=2)



# 1998 to 2007
str(CS98_PH)
str(CS98_PH_QA)

CS98_QA_sel <- CS98_PH %>%
  select(SQUARE_NUM, PLOT_TYPE, REP_NUM, PHF2000) %>%
  mutate(SQUARE_NUM = as.character(SQUARE_NUM)) %>%
  right_join(CS98_PH_QA) %>%
  na.omit()

ggplot(CS98_QA_sel, aes(x = PHF2000, y = PH_DIW_QA)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) 

CS98_QA_sel %>% 
  mutate(PH_DIW_diff = PHF2000 - PH_DIW_QA) %>%
  summarise(diw_mn  = mean(PH_DIW_diff),
            diw_sd  = sd(PH_DIW_diff))
#    diw_mn    diw_sd
# 0.2244382 0.4474434

ph_diff <- CS98_QA_sel$PHF2000 - CS98_QA_sel$PH_DIW_QA
quantile(abs(ph_diff), c(0.95))/1.96
# 0.4920918

ribbon_bounds <- data.frame(x = seq(3.5,9,0.1), 
                            ymax = seq(3.5,9,0.1) + 0.447,
                            ymin = seq(3.5,9,0.1) - 0.447,
                            ymax2 = seq(3.5,9,0.1) + 1.96*0.447,
                            ymin2 = seq(3.5,9,0.1) - 1.96*0.447) 
ggplot() +
  geom_point(data = CS98_QA_sel, aes(x = PHF2000, y = PH_DIW_QA)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin, ymax = ymax), 
              alpha = 0.3, fill = "blue") +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin2, ymax = ymax2), 
              alpha = 0.3, fill = "blue") +
  NULL
# 11 points outside outer ribbon - should have 8.9

ribbon_bounds <- data.frame(x = seq(3.5,9,0.1), 
                            ymax = seq(3.5,9,0.1) + 0.492,
                            ymin = seq(3.5,9,0.1) - 0.492,
                            ymax2 = seq(3.5,9,0.1) + 1.96*0.492,
                            ymin2 = seq(3.5,9,0.1) - 1.96*0.492) 
ggplot() +
  geom_point(data = CS98_QA_sel, aes(x = PHF2000, y = PH_DIW_QA)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin, ymax = ymax), 
              alpha = 0.3, fill = "blue") +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin2, ymax = ymax2), 
              alpha = 0.3, fill = "blue") +
  NULL
# 9 points outside outer ribbon

ribbon_bounds <- data.frame(x = seq(2.5,8.5,0.1), 
                            ymax = seq(2.5,8.5,0.1) + 0.63,
                            ymin = seq(2.5,8.5,0.1) - 0.63,
                            ymax2 = seq(2.5,8.5,0.1) + 1.96*0.63,
                            ymin2 = seq(2.5,8.5,0.1) - 1.96*0.63) 
ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin, ymax = ymax), 
              alpha = 0.3, fill = "blue") +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin2, ymax = ymax2), 
              alpha = 0.3, fill = "blue") +
  # scale_y_continuous(limits = c(2.5,8.5), expand = c(0,0)) +
  scale_x_continuous(limits = c(2.5,8.5), expand = c(0,0)) +
  NULL

MASS::fitdistr(ph_diff, "t")
# m            s            df    
# 0.22816783   0.30264485   3.72142853 
# (0.02716051) (0.02722469) (1.02880468)


# compare visually fits of normal and student T distribution
h <- hist(ph_diff, breaks = 40)
xfit<-seq(min(ph_diff),max(ph_diff),length=40)
yfit<-dnorm(xfit,mean=mean(ph_diff),sd=sd(ph_diff))
yfit <- yfit*diff(h$mids[1:2])*length(ph_diff)
lines(xfit, yfit, col="blue", lwd=2)
yfit<-brms::dstudent_t(xfit,mu = 0.22816783, sigma = 0.30264485, df = 3.72142853)
yfit <- yfit*diff(h$mids[1:2])*length(ph_diff)
lines(xfit, yfit, col="red", lwd=2)



# 2007
CS07_PH_QA_sel <- CS07_PH_QA %>% 
  mutate(REP_ID = paste0(SQUARE_NUM, PLOT_TYPE, REP_NUM)) %>%
  filter(QA_DUP_SAMPLE == "Duplicate Sample") %>%
  select(REP_ID, QA_DUP_SAMPLE, PH_DIW_QA = PH2007_IN_WATER, 
         PH_CACL2_QA = PH2007_IN_CACL2) %>%
  na.omit() %>%
  left_join(select(CS07_PH, REP_ID, PH_DIW = PH2007_IN_WATER,
                   PH_CACL2 = PH2007_IN_CACL2)) %>%
  na.omit()
# 115


ggplot(CS07_PH_QA_sel, aes(x = PH_DIW, y = PH_DIW_QA)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) 
summary(lm(PH_DIW ~ PH_DIW_QA, CS07_PH_QA_sel))

CS07_PH_QA_sel %>% 
  mutate(PH_DIW_diff = PH_DIW - PH_DIW_QA,
         PH_CACL2_diff = PH_CACL2 - PH_CACL2_QA) %>%
  summarise(diw_mn  = mean(PH_DIW_diff),
            diw_sd  = sd(PH_DIW_diff),
            cacl2_mn = mean(PH_CACL2_diff),
            cacl2_sd = sd(PH_CACL2_diff))
#        diw_mn    diw_sd   cacl2_mn  cacl2_sd
#  -0.001478261 0.2527784 0.02156522 0.2273593

ribbon_bounds <- data.frame(x = seq(3.5,9,0.1), 
                            ymax = seq(3.5,9,0.1) + 0.25,
                            ymin = seq(3.5,9,0.1) - 0.25,
                            ymax2 = seq(3.5,9,0.1) + 1.96*0.25,
                            ymin2 = seq(3.5,9,0.1) - 1.96*0.25) 
ggplot() +
  geom_point(data = CS07_PH_QA_sel, aes(x = PH_DIW, y = PH_DIW_QA)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin, ymax = ymax), 
              alpha = 0.3, fill = "blue") +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin2, ymax = ymax2), 
              alpha = 0.3, fill = "blue") +
  NULL
ph07_diff <- CS07_PH_QA_sel$PH_DIW - CS07_PH_QA_sel$PH_DIW_QA
MASS::fitdistr(ph07_diff, "t")
# m            s            df    
# 0.01133412   0.14252233   2.45484216 
# (0.01712093) (0.01895973) (0.70482518)

# compare visually fits of normal and student T distribution
h <- hist(ph07_diff, breaks = 40)
xfit<-seq(min(ph07_diff),max(ph07_diff),length=40)
yfit<-dnorm(xfit,mean=mean(ph07_diff),sd=sd(ph07_diff))
yfit <- yfit*diff(h$mids[1:2])*length(ph07_diff)
lines(xfit, yfit, col="blue", lwd=2)
yfit<-brms::dstudent_t(xfit,mu = 0.01133412, sigma = 0.14252233, df = 2.45484216)
yfit <- yfit*diff(h$mids[1:2])*length(ph07_diff)
lines(xfit, yfit, col="red", lwd=2)


ph07_diff <- CS07_PH_QA_sel$PH_CACL2 - CS07_PH_QA_sel$PH_CACL2_QA
MASS::fitdistr(ph07_diff, "t")
# m            s            df    
# 0.02140830   0.09420440   1.66293369 
# (0.01135743) (0.01451080) (0.39417898)

# 2019
summary(UK19_PH_QA)
UK19_PH_QA %>% as.data.frame() %>%
  mutate(PH_DIW_diff = PH_DIW - PH_DIW_QA,
         PH_CACL2_diff = PH_CACL2 - PH_CACL2_QA) %>%
  summarise(diw_mn  = mean(PH_DIW_diff, na.rm = TRUE),
            diw_sd  = sd(PH_DIW_diff, na.rm = TRUE),
            cacl2_mn = mean(PH_CACL2_diff, na.rm = TRUE),
            cacl2_sd = sd(PH_CACL2_diff, na.rm = TRUE))
#      diw_mn   diw_sd cacl2_mn  cacl2_sd
# -0.02904762 0.165941   -0.015 0.1195322

ribbon_bounds <- data.frame(x = seq(3,8.5,0.1), 
                            ymax = seq(3,8.5,0.1) + 0.16,
                            ymin = seq(3,8.5,0.1) - 0.16,
                            ymax2 = seq(3,8.5,0.1) + 1.96*0.16,
                            ymin2 = seq(3,8.5,0.1) - 1.96*0.16) 
ggplot() +
  geom_point(data = UK19_PH_QA, aes(x = PH_DIW, y = PH_DIW_QA)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin, ymax = ymax), 
              alpha = 0.3, fill = "blue") +
  geom_ribbon(data = ribbon_bounds, aes(x = x, ymin = ymin2, ymax = ymax2), 
              alpha = 0.3, fill = "blue") +
  NULL
# 95% of points within 2 sd of mean

ph19_diff <- na.omit(UK19_PH_QA$PH_DIW - UK19_PH_QA$PH_DIW_QA)
MASS::fitdistr(ph19_diff, "t")
# m             s            df     
# -0.02019703    0.08408127    2.83794645 
# ( 0.01140272) ( 0.01150246) ( 0.94386890)

# compare visually fits of normal and student T distribution
h <- hist(ph19_diff, breaks = 40)
xfit<-seq(min(ph19_diff),max(ph19_diff),length=40)
yfit<-dnorm(xfit,mean=mean(ph19_diff),sd=sd(ph19_diff))
yfit <- yfit*diff(h$mids[1:2])*length(ph19_diff)
lines(xfit, yfit, col="blue", lwd=2)
yfit<-brms::dstudent_t(xfit,mu = -0.02019703, sigma = 0.08408127, df = 2.83794645)
yfit <- yfit*diff(h$mids[1:2])*length(ph19_diff)
lines(xfit, yfit, col="red", lwd=2)


ph19_diff <- na.omit(UK19_PH_QA$PH_CACL2 - UK19_PH_QA$PH_CACL2_QA)
MASS::fitdistr(ph19_diff, "t")
# m              s              df     
# -0.000130451    0.053572120    2.458920628 
# ( 0.007404153) ( 0.007561589) ( 0.748074972)

# compare visually fits of normal and student T distribution
h <- hist(ph19_diff, breaks = 40)
xfit<-seq(min(ph19_diff),max(ph19_diff),length=40)
yfit<-dnorm(xfit,mean=mean(ph19_diff),sd=sd(ph19_diff))
yfit <- yfit*diff(h$mids[1:2])*length(ph19_diff)
lines(xfit, yfit, col="blue", lwd=2)
yfit<-brms::dstudent_t(xfit,mu = -0.000130451, sigma = 0.053572120, df = 2.458920628)
yfit <- yfit*diff(h$mids[1:2])*length(ph19_diff)
lines(xfit, yfit, col="red", lwd=2)


phdiff0719 <- brms::rstudent_t(100000,2.837946,-0.020197,0.084081) -
  brms::rstudent_t(100000,2.454842,0.011334,0.142522)
MASS::fitdistr(phdiff0719, "t")
# m              s              df     
# -0.0308745048    0.1868478133    2.9224355206 
# ( 0.0007269228) ( 0.0007630587) ( 0.0298580039)

phdiff0719 <- brms::rstudent_t(100000,2.458921,-0.000130,0.053572) -
  brms::rstudent_t(100000,1.662934,0.094204,0.021408)
MASS::fitdistr(phdiff0719, "t")
# m               s               df      
# -0.0941307056    0.0651198536    2.3597493212 
# ( 0.0002605570) ( 0.0002760143) ( 0.0208361922)



# Error in rainfall ####
# according to rainfall paper the uncertainty is ~ 0.4%
hist(0.004*ELL_pH$fieldseason_rain)
