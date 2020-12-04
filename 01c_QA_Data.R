# Script for QA analysis for Ellenberg R and soil pH
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())
library(patchwork)

# Ellenberg QA data ####
# Only doing this for Ellenberg R right now
# for small plot calculations can also include Y and U plots (both 4m2)
str(CS_VEG_QA)

CS16_VEG_QA2 <- CS16_VEG_QA %>%
  mutate(REP_ID = substring(Quadrat,2)) %>%
  mutate(PLOT_TYPE = gsub("[0-9]","",REP_ID)) %>%
  filter(PLOT_TYPE %in% c("X","Y","U")) %>%
  select(REP_ID, PLOT_TYPE, BRC_NUMBER = Names, COVER = Cover,
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
  full_join(select(mutate(filter(CS19_SP, 
                                 PLOT_TYPE == "X" &
                                   REP_ID %in% CS19_VEG_QA$REP_ID), 
                          Year = 2019, Surveyor = "CS"),
                   REP_ID, Year, BRC_NUMBER, PLOT_TYPE, COVER = TOTAL_COVER, Surveyor)) %>%
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
  full_join(select(mutate(filter(CS19_SP, 
                                 PLOT_TYPE %in% c("X","XX") & NEST_LEVEL < 2 &
                                   REP_ID %in% CS19_VEG_QA$REP_ID), 
                          Year = 2019, Surveyor = "CS"),
                   REP_ID, Year, BRC_NUMBER, COVER = FIRST_COVER, Surveyor)) %>%
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
  # left_join(select(X_Ell, REP_ID, Year,
  #                  SM_R_W, SM_R_UW,WH_R_W,WH_R_UW) %>%
  #             filter(Year == 2019 & REP_ID %in% filter(X_Ell_QA_inner, Year == 2019)$REP_ID) %>%
  #             mutate(Surveyor = "CS")) %>%
  pivot_longer(contains("_R_"), 
               names_to = c("PlotSize","Ellenberg","Weighting"),
               names_sep = "_", values_to = "score") %>%
  filter(!is.na(score)) %>%
  pivot_wider(names_from = "Surveyor", values_from = "score") %>%
  mutate(Diff = QA - CS)


X_Ell_QA %>%
  mutate(Year = as.factor(Year),
         PlotSize = recode(PlotSize,
                           "SM" = "Small",
                           "WH" = "Full"),
         Weighting = recode(Weighting,
                            "UW" = "Unweighted",
                            "W" = "Weighted")) %>%
  ggplot(aes(x = CS, y = QA, colour = Year)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(Weighting ~ PlotSize) +
  ggthemes::scale_colour_colorblind() +
  coord_fixed() +
  labs(x = "Original survey", y = "QA survey")
ggsave("Ellenberg R QA comparison by plotsize and weighting.png",
       path = "Outputs/Graphs/", width = 15, height = 12, units = "cm")

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
# Year PlotSize Weighting       mu    sd       df
# <dbl> <chr>    <chr>        <dbl> <dbl>    <dbl>
#  1  1990 SM       UW        -0.00521 0.136     2.08
#  2  1990 SM       W         -0.0280  0.197     1.79
#  3  1990 WH       UW        -0.0154  0.124     2.70
#  4  1990 WH       W          0.0424  0.226     4.16
#  5  1998 SM       UW        -0.0593  0.215    17.4 
#  6  1998 SM       W         -0.0269  0.228     3.81
#  7  1998 WH       UW        -0.0980  0.174     8.04
#  8  1998 WH       W         -0.0718  0.239 11587.  
#  9  2007 SM       UW        -0.0174  0.216     3.48
# 10  2007 SM       W          0.00557 0.225     1.79
# 11  2007 WH       UW        -0.0113  0.163     4.05
# 12  2007 WH       W         -0.0393  0.148     2.03
# 13  2019 SM       UW         0.0231  0.248     4.75
# 14  2019 SM       W         -0.0152  0.403     6.39
# 15  2019 WH       UW        -0.0167  0.149     5.68
# 16  2019 WH       W         -0.0406  0.244     2.68

X_Ell_QA_hab <- left_join(X_Ell_QA, BH_IMP)
table(X_Ell_QA_hab$BH_DESC)
table(X_Ell_QA_hab$Management)

ggplot(filter(X_Ell_QA_hab, !is.na(Management)), 
       aes(x = Diff, fill = Management)) + 
  geom_histogram() +
  facet_wrap(~Year + PlotSize + Weighting)

p1 <- X_Ell_QA_hab %>%
  filter(Management == "High") %>%
  ggplot(aes(x = Diff)) + 
  geom_histogram() +
  facet_wrap(~Year + PlotSize + Weighting)
p2 <- X_Ell_QA_hab %>%
  filter(Management == "Low") %>%
  ggplot(aes(x = Diff)) + 
  geom_histogram() +
  facet_wrap(~Year + PlotSize + Weighting)
p1+p2

# seems reasonable to keep the two habitat categories separate

X_Ell_QA_hab <- filter(X_Ell_QA_hab,
                       !is.na(Management))

janitor::tabyl(X_Ell_QA_hab, Year, PlotSize, Management)
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
# Year PlotSize Weighting       mu    sd       df
# <dbl> <chr>    <chr>        <dbl> <dbl>    <dbl>
#   1  1990 SM       UW         0.0407  0.218     5.66
# 2  1990 SM       W         -0.0376  0.264     2.80
# 3  1990 WH       UW        -0.0239  0.118     5.64
# 4  1990 WH       W          0.0450  0.147     2.03
# 5  1998 SM       UW        -0.0662  0.245 17300.  
# 6  1998 SM       W         -0.0475  0.246     4.25
# 7  1998 WH       UW        -0.0947  0.181     8.94
# 8  1998 WH       W         -0.0666  0.249 24044.  
# 9  2007 SM       UW        -0.0185  0.238     4.00
# 10  2007 SM       W          0.0113  0.237     1.98
# 11  2007 WH       UW         0.00441 0.164     3.83
# 12  2007 WH       W         -0.0518  0.161     2.15
# 13  2019 SM       UW         0.0714  0.272     4.50
# 14  2019 SM       W          0.0370  0.444    17.3 
# 15  2019 WH       UW        -0.0247  0.122     5.01
# 16  2019 WH       W         -0.0706  0.325 22478. 

# 1978 (using 1990 data) to 1998
Ell_SMUW_9098 <- brms::rstudent_t(100000,17300,-0.0662,0.245) -
  brms::rstudent_t(100000,5.66,0.0407,0.218)
MASS::fitdistr(Ell_SMUW_9098, "t")
# m              s              df     
# -0.107218165    0.335994946   13.349240068 
# ( 0.001134349) ( 0.001139551) ( 0.454512158)

Ell_SMW_9098 <- brms::rstudent_t(100000,4.25,-0.0475,0.246) -
  brms::rstudent_t(100000,2.8,-0.0376,0.264)
MASS::fitdistr(Ell_SMW_9098, "t")
# m              s              df     
# -0.009222229    0.402272065    3.992830215 
# ( 0.001506537) ( 0.001568434) ( 0.051362460)

Ell_WHUW_9098 <- brms::rstudent_t(100000,8.94,-0.0947,0.181) -
  brms::rstudent_t(100000,5.64,-0.0239,0.118)
MASS::fitdistr(Ell_WHUW_9098, "t")
# m               s               df      
# -0.0713936535    0.2283446128   11.3127031202 
# ( 0.0007785483) ( 0.0008193378) ( 0.3573662196)

Ell_WHW_9098 <- brms::rstudent_t(100000,24044,-0.0666,0.249) -
  brms::rstudent_t(100000,2.03,0.0450,0.147)
MASS::fitdistr(Ell_WHW_9098, "t")
# m              s              df     
# -0.110701734    0.289838786    3.727722331 
# ( 0.001098847) ( 0.001073164) ( 0.042244196)

# 1998 to 2007
Ell_SMUW_9807 <- brms::rstudent_t(100000,4,-0.0185,0.238) - 
  brms::rstudent_t(100000,17300,-0.0662,0.245)
MASS::fitdistr(Ell_SMUW_9807, "t")
# m             s            df     
# 0.047590244   0.345606483   6.755654407 
# (0.001226635) (0.001239306) (0.128356757)

Ell_SMW_9807 <- brms::rstudent_t(100000,1.98,0.0113,0.237) - 
  brms::rstudent_t(100000,4.25,-0.0475,0.246)
MASS::fitdistr(Ell_SMW_9807, "t")
# m             s            df     
# 0.058753608   0.380087881   2.801772237 
# (0.001488855) (0.001542747) (0.027393186)

Ell_WHUW_9807 <- brms::rstudent_t(100000,3.83,0.00441,0.164) - 
  brms::rstudent_t(100000,8.94,-0.0947,0.181)
MASS::fitdistr(Ell_WHUW_9807, "t")
# m              s              df     
# 0.0977445019   0.2589425036   6.6014880982 
# (0.0009206116) (0.0009507893) (0.1268536548)

Ell_WHW_9807 <- brms::rstudent_t(100000,2.15,-0.0518,0.161) - 
  brms::rstudent_t(100000,24044,-0.0666,0.249)
MASS::fitdistr(Ell_WHW_9807, "t")
# m             s            df     
# 0.013821293   0.300038349   3.858779291 
# (0.001131944) (0.001114372) (0.045294892)


# 2007 to 2019
Ell_SMUW_0719 <- brms::rstudent_t(100000,4.5,0.0714,0.272) - 
  brms::rstudent_t(100000,4,-0.0185,0.238) 
MASS::fitdistr(Ell_SMUW_0719, "t")
# m            s            df    
# 0.09095365   0.39628677   5.55613273 
# (0.00143190) (0.00149642) (0.09387955)

Ell_SMW_0719 <- brms::rstudent_t(100000,17.3,0.0370,0.444) - 
  brms::rstudent_t(100000,1.98,0.0113,0.237)
MASS::fitdistr(Ell_SMW_0719, "t")
# m             s            df     
# 0.025829038   0.517221524   3.883124492 
# (0.001950536) (0.001908698) (0.045473849)

Ell_WHUW_0719 <- brms::rstudent_t(100000,5.01,-0.0247,0.122)  - 
  brms::rstudent_t(100000,3.83,0.00441,0.164) 
MASS::fitdistr(Ell_WHUW_0719, "t")
# m               s               df      
# -0.0286482506    0.2224948589    5.2106944507 
# ( 0.0008093376) ( 0.0008399981) ( 0.0827558915)

Ell_WHW_0719 <- brms::rstudent_t(100000,22478,-0.0706,0.325) -
  brms::rstudent_t(100000,2.15,-0.0518,0.161)
MASS::fitdistr(Ell_WHW_0719, "t")
# m              s              df     
# -0.018246730    0.361091203    4.633014585 
# ( 0.001334109) ( 0.001289639) ( 0.061159626)

# Summary stats for s
ELL_SE <- data.frame(
  Time_period = c("7898","9807","0719"),
  ELL_WH_W_SE = c(0.290,0.300,0.361),
  ELL_WH_UW_SE = c(0.228,0.259,0.222),
  ELL_SM_W_SE = c(0.402,0.380,0.517),
  ELL_SM_UW_SE = c(0.336,0.346,0.396)
)



# check normal distribution
X_Ell_QA %>% 
  group_by(Year, PlotSize, Weighting) %>%
  na.omit() %>%
  summarise(mu = MASS::fitdistr(Diff,"normal")$estimate[1],
            sd = MASS::fitdistr(Diff,"normal")$estimate[2])
# Year PlotSize Weighting       mu    sd
# <dbl> <chr>    <chr>        <dbl> <dbl>
#   1  1990 SM       UW        -0.00340 0.240
# 2  1990 SM       W         -0.0764  0.396
# 3  1990 WH       UW         0.00126 0.204
# 4  1990 WH       W          0.0738  0.298
# 5  1998 SM       UW        -0.0572  0.229
# 6  1998 SM       W         -0.0590  0.312
# 7  1998 WH       UW        -0.0907  0.198
# 8  1998 WH       W         -0.0718  0.239
# 9  2007 SM       UW        -0.0384  0.308
# 10  2007 SM       W         -0.0371  0.468
# 11  2007 WH       UW        -0.0251  0.236
# 12  2007 WH       W         -0.0173  0.319
# 13  2019 SM       UW         0.0521  0.320
# 14  2019 SM       W          0.00249 0.480
# 15  2019 WH       UW        -0.00341 0.181
# 16  2019 WH       W          0.00114 0.419

# 1978 (using 1990 data) to 1998
Ell_SMUW_9098 <-rnorm(100000,-0.0572,0.229) -
  rnorm(100000,-0.0034,0.24)
MASS::fitdistr(Ell_SMUW_9098, "normal")
# mean             sd      
# -0.0537842888    0.3331606802 
# ( 0.0010535466) ( 0.0007449699)

Ell_SMW_9098 <- rnorm(100000,-0.059,0.312) -
  rnorm(100000,-0.0764,0.396)
MASS::fitdistr(Ell_SMW_9098, "normal")
# mean           sd     
# 0.016344235   0.502821416 
# (0.001590061) (0.001124343)

Ell_WHUW_9098 <- rnorm(100000,-0.0907,0.198) -
  rnorm(100000,0.00126,0.204)
MASS::fitdistr(Ell_WHUW_9098, "normal")
# mean             sd      
# -0.0933815215    0.2837339246 
# ( 0.0008972455) ( 0.0006344483)

Ell_WHW_9098 <- rnorm(100000,-0.0718,0.239) -
  rnorm(100000,0.0738,0.298)
MASS::fitdistr(Ell_WHW_9098, "normal")
# mean             sd      
# -0.1452812734    0.3823796764 
# ( 0.0012091907) ( 0.0008550269)

# 1998 to 2007
Ell_SMUW_9807 <- rnorm(100000,-0.0384,0.308) - 
  rnorm(100000,-0.0572,0.229)
MASS::fitdistr(Ell_SMUW_9807, "normal")
# mean            sd     
# 0.0188350164   0.3848821075 
# (0.0012171041) (0.0008606226)

Ell_SMW_9807 <- rnorm(100000,-0.0371,0.468) - 
  rnorm(100000,-0.0590,0.312)
MASS::fitdistr(Ell_SMW_9807, "normal")
# mean           sd     
# 0.022972157   0.562380581 
# (0.001778404) (0.001257521)

Ell_WHUW_9807 <- rnorm(100000,-0.0251,0.236) - 
  rnorm(100000,-0.0907,0.198)
MASS::fitdistr(Ell_WHUW_9807, "normal")
# mean            sd     
# 0.0660062989   0.3077434036 
# (0.0009731701) (0.0006881352)

Ell_WHW_9807 <- rnorm(100000,-0.0173,0.319) - 
  rnorm(100000,-0.0718,0.239)
MASS::fitdistr(Ell_WHW_9807, "normal")
# mean            sd     
# 0.0545133358   0.3973044908 
# (0.0012563871) (0.0008883998)


# 2007 to 2019
Ell_SMUW_0719 <- rnorm(100000,0.0521,0.320) - 
  rnorm(100000,-0.0384,0.308) 
MASS::fitdistr(Ell_SMUW_0719, "normal")
# mean            sd     
# 0.0901869514   0.4431574876 
# (0.0014013870) (0.0009909303)

Ell_SMW_0719 <- rnorm(100000,0.00249,0.480) - 
  rnorm(100000,-0.0371,0.468)
MASS::fitdistr(Ell_SMW_0719, "normal")
# mean           sd     
# 0.040685935   0.672016834 
# (0.002125104) (0.001502675)

Ell_WHUW_0719 <- rnorm(100000,-0.00341,0.181)  - 
  rnorm(100000,-0.0251,0.236) 
MASS::fitdistr(Ell_WHUW_0719, "normal")
# mean            sd     
# 0.0222938300   0.2972070206 
# (0.0009398511) (0.0006645751)

Ell_WHW_0719 <- rnorm(100000,0.00114,0.419) -
  rnorm(100000,-0.0173,0.319)
MASS::fitdistr(Ell_WHW_0719, "normal")
# mean           sd     
# 0.018484580   0.525938009 
# (0.001663162) (0.001176033)

# Summary stats for s
ELL_SE <- data.frame(
  Time_period = c("7898","9807","0719"),
  ELL_WH_W_SE_NORM = c(0.382,0.397,0.526),
  ELL_WH_UW_SE_NORM = c(0.284,0.308,0.297),
  ELL_SM_W_SE_NORM = c(0.503,0.562,0.672),
  ELL_SM_UW_SE_NORM = c(0.333,0.385,0.443)
)


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

MASS::fitdistr(ph78_diff, "normal")
# mean           sd     
# -0.27727011    0.46097491 
# ( 0.03494642) ( 0.02471085)

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

MASS::fitdistr(ph_diff, "normal")
# mean          sd    
# 0.22443820   0.44618479 
# (0.03344296) (0.02364775)

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
MASS::fitdistr(ph07_diff, "normal")
# mean            sd     
# -0.001478261    0.251676962 
# ( 0.023468998) ( 0.016595087)

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
MASS::fitdistr(ph07_diff, "normal")
# mean          sd    
# 0.02156522   0.22636865 
# (0.02110899) (0.01492631)

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
MASS::fitdistr(ph19_diff, "normal")
# mean           sd     
# -0.02904762    0.16495035 
# ( 0.01799756) ( 0.01272620)

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
MASS::fitdistr(ph19_diff, "normal")
# mean            sd     
# -0.015000000    0.118818589 
# ( 0.012964171) ( 0.009167053)

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

phdiff0719 <- rnorm(100000,-0.02904762,0.16495035) -
  rnorm(100000,-0.001478261,0.251676962)
MASS::fitdistr(phdiff0719, "normal")
# mean             sd      
# -0.0273476998    0.3027112396 
# ( 0.0009572570) ( 0.0006768829)

phdiff0719 <- rnorm(100000,-0.015,0.118818589) -
  rnorm(100000,0.02156522,0.22636865)
MASS::fitdistr(phdiff0719, "normal")
# mean             sd      
# -0.0374323454    0.2565189188 
# ( 0.0008111840) ( 0.0005735937)