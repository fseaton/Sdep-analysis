# Script for various data manipulations
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(sf)

# Repeat plot identification ####
# first set up pH REP_IDs
UK19_PH <- UK19_PH %>%
  mutate(REP_ID = paste(SQUARE_NUM,REP_NUM, sep = "X")) %>%
  mutate(REP_ID = recode(REP_ID,
                         "546X3" = "546X6",
                         "546X4" = "546X7",
                         "776X5" = "776X6",
                         "828X8" = "828X7",
                         "861X4" = "861X6",
                         "878X2" = "878X6",
                         "983X2" = "983X6",
                         "1056X3" = "1056X6"))


CS78_PH$REP_ID <- paste(CS78_PH$SQUARE_NUM,CS78_PH$REP_NUM, sep = "X")
CS98_PH$REP_ID <- paste(CS98_PH$SQUARE_NUM,CS98_PH$REP_NUM, sep = "X")
CS07_PH$REP_ID <- paste(CS07_PH$SQUARE_NUM,CS07_PH$REP_NUM, sep = "X")
CS16_PH$REP_ID <- paste(CS16_PH$SQUARE_NUM,CS16_PH$REP_NUM, sep = "X")

# now correct repeat plots for 2019 data and missing soils plots
CS_REP_ID <- filter(CS_REPEAT_PLOTS, AMALG_PTYPE == "X") %>%
  mutate(Y19 = ifelse(!is.na(Y07), Y07,
                      ifelse(!is.na(Y9899), Y9899,
                             ifelse(!is.na(Y90), Y90,
                                    Y78))),
         Y78 = ifelse(is.na(Y78) & !Y78 %in% CS78_PH$REP_ID & Y90 %in% CS78_PH$REP_ID, 
                      Y90, 
                      ifelse(!Y78 %in% CS78_PH$REP_ID & Y9899 %in% CS78_PH$REP_ID, Y9899, Y78)),
         Y9899 = ifelse(is.na(Y9899) & !Y9899 %in% CS98_PH$REP_ID & Y90 %in% CS98_PH$REP_ID, 
                              Y90, Y9899)) %>%
  mutate(Y07 = ifelse(is.na(Y07) & !Y07 %in% CS07_PH$REP_ID & Y19 %in% CS07_PH$REP_ID,
                      Y19, Y07))


# AVC data manipulation ####
hab07 <- select(CS07_IBD, REP_ID = REP_ID07, AVC07) %>%
  unique()
hab98 <- select(CS98_IBD, REP_ID = REP_ID98, AVC98) %>%
  unique()
hab90 <- select(CS90_IBD, REP_ID = REP_ID90, AVC90) %>%
  unique()
hab78 <- select(CS78_IBD, REP_ID = REP_ID78, AVC78) %>%
  unique()

# create combined AVC variable, if 07 has AVC use that otherwise use 98 then 78.
# There are only 3 sites with no AVC data and I can't see how to get theirs as
# they don't appear in 2016/19.
hab <- full_join(hab07, hab98) %>% 
  full_join(hab90) %>% full_join(hab78) %>%
  mutate_if(is.factor, as.character) %>%
  mutate(AVC = ifelse(!is.na(AVC07), AVC07,
                      ifelse(!is.na(AVC98), AVC98,
                             ifelse(!is.na(AVC90), AVC90,
                                    ifelse(!is.na(AVC78), AVC78, NA))))) %>%
  mutate(AVC_desc = recode(AVC,
                           `1` = "Crops/Weeds",
                           `2` = "Tall herb/grass",
                           `3` = "Fertile grassland",
                           `4` = "Infertile grassland",
                           `5` = "Lowland wooded",
                           `6` = "Upland wooded",
                           `7` = "Moorland grass/mosaic",
                           `8` = "Heath/bog"))

# Broad habitat data ####
BH_comb <- do.call(rbind,
                   list(mutate(select(VEGETATION_PLOTS_20161819, 
                                      REP_ID, BH = BH_PLOT),
                               Year = 2019) %>%
                          left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y19)) %>%
                          mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), 
                                                 REPEAT_PLOT_ID, REP_ID)) %>%
                          select(-REPEAT_PLOT_ID),
                        mutate(select(CS07_IBD, 
                                      REP_ID = REP_ID07, BH = BH07),
                               Year = 2007) %>%
                          left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y07)) %>%
                          mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), 
                                                 REPEAT_PLOT_ID, REP_ID)) %>%
                          select(-REPEAT_PLOT_ID),
                        mutate(select(CS98_IBD, 
                                      REP_ID = REP_ID98, BH = BH98),
                               Year = 1998) %>%
                          left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y9899)) %>%
                          mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), 
                                                 REPEAT_PLOT_ID, REP_ID)) %>%
                          select(-REPEAT_PLOT_ID),
                        mutate(select(CS90_IBD, 
                                      REP_ID = REP_ID90, BH = BH90),
                               Year = 1990)%>%
                          left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y90)) %>%
                          mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), 
                                                 REPEAT_PLOT_ID, REP_ID)) %>%
                          select(-REPEAT_PLOT_ID),
                        mutate(select(CS78_IBD, 
                                      REP_ID = REP_ID78, BH = BH78),
                               Year = 1978)%>%
                          left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y78)) %>%
                          mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), 
                                                 REPEAT_PLOT_ID, REP_ID)) %>%
                          select(-REPEAT_PLOT_ID))) %>%
  left_join(unique(select(BHPH_NAMES, BH = BH_CODE,
                          BH_DESC = BROAD_HABITAT)))

BH_comb_nodupes <- BH_comb %>%
  group_by(REP_ID, Year) %>%
  summarise(BH = min(BH), .groups = "drop") %>%
  left_join(unique(select(BHPH_NAMES, BH = BH_CODE,
                          BH_DESC = BROAD_HABITAT)))
  


# Ellenberg scores ####
str(CS19_SP)
table(CS19_SP$PLOT_TYPE)
unique(CS19_SP[CS19_SP$PLOT_TYPE=="XX","REP_ID"])
table(CS19_SP[CS19_SP$PLOT_TYPE=="X","PLOTYEAR"])

str(SPECIES_LIB_TRAITS)
filter(SPECIES_LIB_CODES, COLUMN_NAME == "GROWTH_FORM")

CS18_ELL <- filter(CS19_SP, PLOT_TYPE %in% c("X","XX")) %>%
  mutate(REP_ID = paste0(SQUARE,PLOT_TYPE,PLOT_NUMBER)) %>%
  mutate(REP_ID = gsub("XX","X",REP_ID)) %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   starts_with("EBER"),
                   GROWTH_FORM)) %>%
  filter(GROWTH_FORM %in% c("f","fe","g","m","s","ss","w")) %>% # filter to vascular plants
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% # set 0 values to NA
  group_by(REP_ID) %>%
  summarise(across(starts_with("EBER"), mean, na.rm = TRUE,
                   .names = "{col}18")) %>%
  rename_with(~gsub("EBERG","",.x))
summary(CS18_ELL)
test <- full_join(CS18_ELL, GM16_IBD, by = c("REP_ID" = "REP_ID16"))
plot(N18 ~ N16, test);abline(0,1)

CS98_ELL <- CS98_SP %>%
  select(REP_ID, BRC_NUMBER, TOTAL_COVER) %>%
  unique() %>%
  filter(TOTAL_COVER > 0) %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   starts_with("EBER"),
                   GROWTH_FORM)) %>%
  # filter(GROWTH_FORM %in% c("f","fe","g","m","s","ss","w")) %>% # filter to vascular plants
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% # set 0 values to NA
  group_by(REP_ID) %>%
  summarise(across(starts_with("EBER"), function(x) sum(x, na.rm=TRUE)/length(na.omit(EBERGN)),
                   .names = "{col}98_new")) %>%
  rename_with(~gsub("EBERG","",.x))
test <- full_join(CS98_ELL, CS98_IBD, by = c("REP_ID" = "REP_ID98"))
#par(mfrow=c(2,2))
plot(R98_new ~ R98, test);abline(0,1)
plot(N98_new ~ N98, test);abline(0,1)
plot(W98_new ~ F98, test);abline(0,1)
plot(L98_new ~ L98, test);abline(0,1)
par(mfrow=c(1,1))
summary(CS18_ELL)

# Ellenberg for inner 2x2m square
X_Ell_inner <- CS07_SP %>%
  left_join(select(CS07_PLOTS, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
  filter(PLOT_TYPE == "X" & NEST_LEVEL < 2)  %>%
  mutate(Year = 2007) %>%
  select(REP_ID, Year, BRC_NUMBER, FIRST_COVER, TOTAL_COVER) %>%
  full_join(select(mutate(filter(CS19_SP, 
                                 PLOT_TYPE %in% c("X","XX") & NEST_LEVEL < 2), Year = 2019),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER, TOTAL_COVER)) %>%
  full_join(select(mutate(filter(CS98_SP, grepl("X", REP_ID) & NEST_LEVEL < 2), Year = 1998),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER, TOTAL_COVER)) %>%
  full_join(select(mutate(filter(CS90_SP, grepl("X", REP_ID) & QUADRAT_NEST < 2), Year = 1990),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER = C1, TOTAL_COVER)) %>%
  full_join(select(mutate(filter(CS78_SP, grepl("X", REP_ID) & QUADRAT_NEST < 2), Year = 1978),
                   REP_ID, Year, BRC_NUMBER, TOTAL_COVER = COVER)) %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   starts_with("EBER"))) %>%
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% # set 0 values to NA
  group_by(Year, REP_ID) %>%
  summarise(across(starts_with("EBER"), mean, na.rm = TRUE)) %>%
  rename_with(~gsub("EBERG","SM_",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA))

X_Ell_whole <- CS07_SP %>%
  left_join(select(CS07_PLOTS, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
  filter(PLOT_TYPE == "X")  %>%
  mutate(Year = 2007) %>%
  select(REP_ID, Year, BRC_NUMBER, TOTAL_COVER) %>% unique() %>%
  full_join(unique(select(mutate(filter(CS19_SP, 
                                        PLOT_TYPE == "X"), Year = 2019),
                          REP_ID, Year, BRC_NUMBER, TOTAL_COVER))) %>%
  full_join(unique(select(mutate(filter(CS98_SP, grepl("X", REP_ID)), Year = 1998),
                          REP_ID, Year, BRC_NUMBER, TOTAL_COVER))) %>%
  full_join(unique(select(mutate(filter(CS90_SP, grepl("X", REP_ID)), Year = 1990),
                          REP_ID, Year, BRC_NUMBER, TOTAL_COVER))) %>%
  full_join(unique(select(mutate(filter(CS78_SP, grepl("X", REP_ID)), Year = 1978),
                          REP_ID, Year, BRC_NUMBER, TOTAL_COVER = COVER))) %>%
  left_join(unique(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                          starts_with("EBER")))) %>%
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% # set 0 values to NA
  group_by(Year, REP_ID) %>%
  summarise(across(starts_with("EBER"), mean, na.rm = TRUE)) %>%
  rename_with(~gsub("EBERG","WH_",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA))

X_Ell_comp <- full_join(X_Ell_inner, X_Ell_whole) %>%
  left_join(hab) %>%
  mutate(R_diff = SM_R - WH_R,
         N_diff = SM_N - WH_N,
         W_diff = SM_W - WH_W,
         L_diff = SM_L - WH_L)

p1 <- ggplot(X_Ell_comp, aes(x = WH_R, y = SM_R)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc) +
  labs(x = bquote("Ellenberg R 400m"^2~"plot"),
       y = bquote("Ellenberg R 4m"^2~"plot"))

ggplot(X_Ell_comp, aes(x = WH_N, y = SM_N)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc)

ggplot(X_Ell_comp, aes(x = R_diff)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~AVC_desc)

ggplot(X_Ell_comp, aes(x = WH_R, y = SM_R)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(Year~AVC_desc) +
  theme_bw()

ggplot(X_Ell_comp, aes(x = N_diff)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~AVC_desc)

X_Ell_comp %>%
  select(Year, REP_ID, ends_with("diff"), AVC_desc) %>%
  pivot_longer(ends_with("diff"), names_to = "Ellenberg") %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-2.5,2.5)) +
  facet_grid(Ellenberg~AVC_desc)

# weighted mean Ellenberg 
X_wEll_inner <- CS07_SP %>%
  left_join(select(CS07_PLOTS, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
  filter(PLOT_TYPE == "X" & NEST_LEVEL < 2)  %>%
  mutate(Year = 2007) %>%
  select(REP_ID, Year, BRC_NUMBER, FIRST_COVER, TOTAL_COVER) %>%
  full_join(select(mutate(filter(CS19_SP, 
                                 PLOT_TYPE %in% c("X","XX") & NEST_LEVEL < 2), Year = 2019),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER, TOTAL_COVER)) %>%
  full_join(select(mutate(filter(CS98_SP, grepl("X", REP_ID) & NEST_LEVEL < 2), Year = 1998),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER, TOTAL_COVER)) %>%
  full_join(select(mutate(filter(CS90_SP, grepl("X", REP_ID) & QUADRAT_NEST < 2), Year = 1990),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER = C1, TOTAL_COVER)) %>%
  full_join(select(mutate(filter(CS78_SP, grepl("X", REP_ID)), Year = 1978),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER = COVER)) %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   starts_with("EBER"))) %>%
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% 
  group_by(Year, REP_ID) %>%
  summarise(across(starts_with("EBER"), weighted.mean, w = FIRST_COVER, na.rm = TRUE)) %>%
  rename_with(~gsub("EBERG","SM_",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA))

X_wEll_whole <- CS07_SP %>%
  left_join(select(CS07_PLOTS, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
  filter(PLOT_TYPE == "X")  %>%
  mutate(Year = 2007) %>%
  select(REP_ID, Year, BRC_NUMBER, TOTAL_COVER) %>% unique() %>%
  full_join(unique(select(mutate(filter(CS19_SP, 
                                        PLOT_TYPE == "X"), Year = 2019),
                          REP_ID, Year, BRC_NUMBER, TOTAL_COVER))) %>%
  full_join(unique(select(mutate(filter(CS98_SP, grepl("X", REP_ID)), Year = 1998),
                          REP_ID, Year, BRC_NUMBER, TOTAL_COVER))) %>%
  full_join(unique(select(mutate(filter(CS90_SP, grepl("X", REP_ID)), Year = 1990),
                          REP_ID, Year, BRC_NUMBER, TOTAL_COVER))) %>%
  full_join(unique(select(mutate(filter(CS78_SP, grepl("X", REP_ID)), Year = 1978),
                          REP_ID, Year, BRC_NUMBER, TOTAL_COVER = COVER))) %>%
  left_join(unique(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                          starts_with("EBER")))) %>%
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% 
  group_by(Year, REP_ID) %>%
  summarise(across(starts_with("EBER"), weighted.mean, w = TOTAL_COVER, na.rm = TRUE)) %>%
  rename_with(~gsub("EBERG","WH_",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA))

X_wEll_comp <- full_join(X_wEll_inner, X_wEll_whole) %>%
  left_join(hab) %>%
  mutate(R_diff = SM_R - WH_R,
         N_diff = SM_N - WH_N,
         W_diff = SM_W - WH_W,
         L_diff = SM_L - WH_L)

p2 <- ggplot(filter(X_wEll_comp, Year != 1978), aes(x = WH_R, y = SM_R)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc) +
  labs(x = bquote("Ellenberg R 400m"^2~"plot"),
       y = bquote("Ellenberg R 4m"^2~"plot"))

ggplot(filter(X_wEll_comp, Year != 1978), aes(x = R_diff)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~AVC_desc)

X_wEll_comp %>%
  select(Year, REP_ID, ends_with("diff"), AVC_desc) %>%
  pivot_longer(ends_with("diff"), names_to = "Ellenberg") %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-2.5,2.5)) +
  facet_grid(Ellenberg~AVC_desc)


t.test(X_wEll_comp$SM_R,X_wEll_comp$WH_R)
# t = -0.47853, df = 18816, p-value = 0.6323
t.test(X_Ell_comp$SM_R,X_Ell_comp$WH_R)
# t = -3.3001, df = 19015, p-value = 0.0009682

x <- na.omit(unique(X_wEll_comp$AVC_desc))
for(i in 1:length(x)){
  
  dat <- filter(X_wEll_comp, AVC_desc == x[i])
  print(paste(x[i],"p =",round(t.test(dat$SM_R,dat$WH_R)$p.value,5)))
}
# [1] "Tall herb/grass p = 0.85153"
# [1] "Crops/Weeds p = 0.04266"
# [1] "Fertile grassland p = 0.92217"
# [1] "Heath/bog p = 2e-05"
# [1] "Moorland grass/mosaic p = 0.12709"
# [1] "Upland wooded p = 0.74556"
# [1] "Infertile grassland p = 0.92811"
# [1] "Lowland wooded p = 0.39651"

x <- na.omit(unique(X_Ell_comp$AVC_desc))
for(i in 1:length(x)){
  dat <- filter(X_Ell_comp, AVC_desc == x[i])
  print(paste(x[i],"p =",round(t.test(dat$SM_R,dat$WH_R)$p.value,5)))
}
# [1] "Tall herb/grass p = 0.91431"
# [1] "Crops/Weeds p = 0.06944"
# [1] "Fertile grassland p = 0.05292"
# [1] "Heath/bog p = 0"
# [1] "Moorland grass/mosaic p = 0"
# [1] "Upland wooded p = 0.00329"
# [1] "Infertile grassland p = 0.06137"
# [1] "Lowland wooded p = 0.96831"

p1 + ggtitle("Unweighted") + p2 + ggtitle("Cover weighted")
ggsave("Ellenberg R plot size comparison.png", path = "Outputs/Graphs/",
       width = 24, height = 12, units = "cm")

p1 <- ggplot(filter(X_Ell_comp, Year != 1978), aes(x = WH_N, y = SM_N)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc) +
  labs(x = bquote("Ellenberg N 400m"^2~"plot"),
       y = bquote("Ellenberg N 4m"^2~"plot")) +
  ggtitle("Unweighted")
p2 <- ggplot(filter(X_wEll_comp, Year != 1978), aes(x = WH_N, y = SM_N)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc) +
  labs(x = bquote("Ellenberg N 400m"^2~"plot"),
       y = bquote("Ellenberg N 4m"^2~"plot")) +
  ggtitle("Cover weighted")
p1 + p2
ggsave("Ellenberg N plot size comparison.png", path = "Outputs/Graphs/",
       width = 24, height = 12, units = "cm")

t.test(X_wEll_comp$SM_N,X_wEll_comp$WH_N)
# t = -0.12149, df = 18823, p-value = 0.9033
t.test(X_Ell_comp$SM_N,X_Ell_comp$WH_N)
# t = -1.621, df = 19043, p-value = 0.105

# correlations
x <- na.omit(unique(X_wEll_comp$AVC_desc))
for(i in 1:length(x)){
  
  dat <- filter(X_wEll_comp, AVC_desc == x[i])
  print(paste(x[i],"p =",round(t.test(dat$SM_N,dat$WH_N)$p.value,5)))
}
# [1] "Tall herb/grass p = 0.8639"
# [1] "Crops/Weeds p = 0.04501"
# [1] "Fertile grassland p = 0.83626"
# [1] "Heath/bog p = 0.03957"
# [1] "Moorland grass/mosaic p = 0.07602"
# [1] "Upland wooded p = 0.63365"
# [1] "Infertile grassland p = 0.38431"
# [1] "Lowland wooded p = 0.31157"

x <- na.omit(unique(X_Ell_comp$AVC_desc))
for(i in 1:length(x)){
  dat <- filter(X_Ell_comp, AVC_desc == x[i])
  print(paste(x[i],"p =",round(t.test(dat$SM_N,dat$WH_N)$p.value,5)))
}
# [1] "Tall herb/grass p = 0.91795"
# [1] "Crops/Weeds p = 0.05538"
# [1] "Fertile grassland p = 0.38703"
# [1] "Heath/bog p = 0"
# [1] "Moorland grass/mosaic p = 0"
# [1] "Upland wooded p = 0.01399"
# [1] "Infertile grassland p = 0.41287"
# [1] "Lowland wooded p = 0.62069"

# all data
X_Ell <- full_join(
  rename_with(X_wEll_inner, ~ paste0(.x, "_W"), starts_with("SM")),
  rename_with(X_wEll_whole, ~ paste0(.x, "_W"), starts_with("WH"))
) %>%
  full_join(
    rename_with(X_Ell_inner, ~ paste0(.x, "_UW"), starts_with("SM"))
  ) %>%
  full_join(
    rename_with(X_Ell_whole, ~ paste0(.x, "_UW"), starts_with("WH"))
  )

# Ellenberg QA data ####
# Only doing this for Ellenberg R right now
str(CSVEG_QA)

X_Ell_whole_QA <- CSVEG_QA %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   EBERGR)) %>%
  mutate(across(EBERGR, na_if, y = 0)) %>% # set 0 values to NA
  group_by(Year, Surveyor, REP_ID) %>%
  summarise(across(EBERGR, mean, na.rm = TRUE)) %>%
  rename_with(~gsub("EBERG","WH_",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA)) %>%
  pivot_wider(names_from = "Surveyor", values_from = "WH_R") %>%
  mutate(WH_UW_diff = CS - QA)

ggplot(X_Ell_whole_QA, aes(x = CS, y = QA, colour = Year)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

EllR_diff <- na.omit(X_Ell_whole_QA$WH_UW_diff)

X_Ell_whole_QA %>% 
  mutate(Year = ifelse(Year == 2016, 2019, Year)) %>%
  group_by(Year) %>%
  na.omit() %>%
  summarise(mu = MASS::fitdistr(WH_UW_diff,"t")$estimate[1],
            sd = MASS::fitdistr(WH_UW_diff,"t")$estimate[2],
            df = MASS::fitdistr(WH_UW_diff,"t")$estimate[3])
# Year      mu    sd    df
# <dbl>   <dbl> <dbl> <dbl>
# 1  1990 0.00775 0.121  2.58
# 2  1998 0.103   0.169  6.95
# 3  2007 0.0111  0.164  4.08
# 4  2019 0.0169  0.150  5.25

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


# pH data ####
# Quick fix for UK19_PH values that aren't matching to the veg plots
UK19_PH <- UK19_PH %>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y19)) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID)) %>%
  select(-REPEAT_PLOT_ID)

CS78_PH <- CS78_PH %>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y78)) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID))

CS98_PH <- CS98_PH %>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y9899)) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID))

CS07_PH <- CS07_PH %>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y07)) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID))

CS16_PH <- CS16_PH %>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y19)) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID))

# wide format pH
PH <- full_join(select(CS78_PH, REP_ID, PH_1978 = PH1978),
                select(CS98_PH, REP_ID, PH_1998 = PHF2000)) %>%
  full_join(select(CS07_PH, REP_ID, PH_2007 = PH2007_IN_WATER, PHC_2007 = PH2007_IN_CACL2)) %>%
  full_join(select(CS16_PH, REP_ID, PH_2016 = PH_DIW, PHC_2016 = PH_CACL2)) %>%
  full_join(select(UK19_PH, REP_ID, PH_2019 = PH_DIW, PHC_2019 = PH_CACL))
str(PH)
summary(PH)
mice::md.pattern(PH)

# long format pH
PH_long <- pivot_longer(PH, starts_with("PH"),
                        values_to = "pH",
                        values_drop_na = TRUE,
                        names_to = c("variable","Year"),
                        names_sep = "_",
                        names_transform = list(Year = as.integer)) %>%
  mutate(variable = ifelse(variable == "PH", "pH","pH_CaCl2")) %>%
  pivot_wider(names_from = variable, values_from = pH)
PH_long

# pH differences
# calculate differences between survey years 
PH2 <- PH_long %>%
  left_join(CS_REP_ID) %>%
  select(REP_ID = REPEAT_PLOT_ID)
  mutate(diff7898 = PH_1998 - PH_1978,
         diff7807 = PH_2007 - PH_1978,
         diff7816 = PH_2016 - PH_1978,
         diff7819 = PH_2019 - PH_1978,
         diff9807 = PH_2007 - PH_1998,
         diff9816 = PH_2016 - PH_1998,
         diff9819 = PH_2019 - PH_1998,
         diff0716 = PH_2016 - PH_2007,
         diff0719 = PH_2019 - PH_2007) %>%
  mutate(diff0718 = ifelse(!is.na(diff0719), diff0719,
                           ifelse(!is.na(diff0716), diff0716, NA)),
         diff7818 = ifelse(!is.na(diff7819), diff7819,
                           ifelse(!is.na(diff7816), diff7816, NA)),
         diff9818 = ifelse(!is.na(diff9819), diff9819,
                           ifelse(!is.na(diff9816), diff9816, NA)))
summary(PH)

PH_diff_long <- PH %>%
  select(REP_ID, starts_with("diff")) %>%
  pivot_longer(starts_with("diff"),
               values_to = "pH",
               values_drop_na = TRUE) %>%
  mutate(name = as.factor(name)) %>%
  mutate(name = forcats::fct_inorder(name))

# select only most recent change and convert into wide format for plotting
PH_Diff_wide <- select(PH, REP_ID, diff0718) %>%
  na.omit() %>%
  left_join(select(CS07_PLOTS, REP_ID, POINT_X, POINT_Y))


# LOI ####
LOI <- full_join(select(CS78_PH, REP_ID, LOI_1978 = LOI1978),
                select(CS98_PH, REP_ID, LOI_1998 = LOI2000)) %>%
  full_join(select(CS07_PH, REP_ID, LOI_2007 = LOI2007)) %>%
  full_join(select(UK19_PH, REP_ID, LOI_2019 = LOI)) %>%
  mutate(LOI_2019 = 100*LOI_2019)
str(LOI)
summary(LOI)

# long format pH
LOI_long <- pivot_longer(LOI, starts_with("LOI"),
                        values_to = "LOI",
                        values_drop_na = TRUE,
                        names_to = "Year",
                        names_prefix = "LOI_",
                        names_transform = list(Year = as.integer)) 
LOI_long


# Soil moisture ####
UK19_WET <- UK19_WET %>%
  mutate(REP_ID = paste0(`Square number...1`,`X plot...2`)) %>%
  select(-starts_with("Square")) %>%
  select(-starts_with("X "))

MOISTURE <- CS_tier4 %>%
  mutate(REP_ID = ifelse(!is.na(REP_ID07), REP_ID07,
                         ifelse(!is.na(REP_ID98), REP_ID98,
                                ifelse(!is.na(REP_ID78), REP_ID78, NA)))) %>%
  select(REP_ID, MOISTURE_CONTENT_07, MOISTURE_CONTENT_98) %>%
  full_join(select(UK19_WET, REP_ID, MOISTURE_CONTENT_19 = `g water/wet weight of soil`)) %>%
  # mutate(VWC_19 = 100*VWC_19) %>%
  pivot_longer(starts_with("MOISTURE"), names_to = "variable", 
               values_to = "Moisture") %>%
  mutate(Year = ifelse(variable == "MOISTURE_CONTENT_07", 2007,
                       ifelse(variable == "MOISTURE_CONTENT_98", 1998,
                              ifelse(variable == "MOISTURE_CONTENT_19", 2019, NA)))) %>%
  select(-variable) %>% na.omit()


# Spatial data manipulation ####
# rainfall - monthly totals ####
rain07 <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_200701-200712.nc"
rain <- raster::brick(rain07)
rain[rain > 9e20] <- NA

cs_loc07 <- plot_locations %>%
  filter(YEAR == "y07") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain <- raster::extract(rain, cs_loc07)
rownames(cs_loc_rain) <- cs_loc07$REP_ID
# for reps that don't fall in the raster get average of all values within 2000m
cs_loc07_na <- cs_loc07 %>%
  filter(REP_ID %in% rownames(cs_loc_rain)[is.na(cs_loc_rain[,1])])
cs_loc_rain_na <- raster::extract(rain, cs_loc07_na,
                                  fun = mean, buffer = 2000)
rownames(cs_loc_rain_na) <- cs_loc07_na$REP_ID
cs_loc_rain <- rbind(na.omit(cs_loc_rain),
                     cs_loc_rain_na)
colnames(cs_loc_rain) <- paste("X",c(1:12), "2007", sep = "_")


cs_loc_rain07_long <- cs_loc_rain %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  pivot_longer(starts_with("X"), names_to = c("Month", "Year"),
               names_prefix = "X_", names_sep = "_",
               values_to = "rainfall") %>%
  mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "[A-Z]"),"[", 1)),
         rainfall = ifelse(rainfall < 9e20, rainfall, NA))

str(plot_dates)
sample_date <- plot_dates %>%
  select(SERIES_NUM, starts_with("MID_")) %>%
  mutate(across(starts_with("MID_"), lubridate::ymd)) %>%
  mutate(across(starts_with("MID_"), lubridate::month))

cs_loc_rain07_long <- cs_loc_rain07_long %>%
  left_join(select(sample_date, SERIES_NUM, DATE = MID__DATE07)) %>%
  mutate(DATE = ifelse(!is.na(DATE), DATE, 
                       round(mean(DATE, na.rm = TRUE))+1),
         Month = as.numeric(Month)) %>%
  mutate(field_season = ifelse(Month <= DATE & Month > DATE - 4, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_rainfall = mean(rainfall),
            sum_rainfall = sum(rainfall)) 


# 1998
rain98 <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_199801-199812.nc"
rain <- raster::brick(rain98)
rain[rain > 9e20] <- NA

# check to see plots and rasters align
cs_loc98_wgs <- plot_locations %>%
  filter(YEAR == "y9899") %>% 
  mutate(EAST = ifelse(!is.na(E_6_FIG_100M),E_6_FIG_100M, E_4_FIG_1KM),
         NORTH = ifelse(!is.na(N_6_FIG_100M),N_6_FIG_100M, N_4_FIG_1KM)) %>% 
  select(SERIES_NUM, EAST, NORTH) %>%
  group_by(SERIES_NUM) %>%
  summarise(EAST = mean(EAST), NORTH = mean(NORTH)) %>%
  ungroup() %>%
  mutate(EAST= round(EAST, -3), NORTH = round(NORTH, -3)) %>%
  na.omit() %>%
  st_as_sf(coords = c("EAST","NORTH"), crs = 27700) %>%
  st_transform(4326)
leaflet() %>% addTiles() %>% 
  addRasterImage(rain_june, opacity = 0.5) %>% 
  addMarkers(data=cs_loc98_wgs, 
             label = cs_loc98_wgs$SERIES_NUM)

cs_loc98 <- plot_locations %>%
  filter(YEAR == "y9899") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain98 <- raster::extract(rain, cs_loc98)
rownames(cs_loc_rain98) <- cs_loc98$REP_ID
# for reps that don't fall in the raster get average of all values within 2000m
cs_loc98_na <- cs_loc98 %>%
  filter(REP_ID %in% rownames(cs_loc_rain)[is.na(cs_loc_rain[,1])])
cs_loc_rain_na <- raster::extract(rain, cs_loc98_na,
                                  fun = mean, buffer = 2000)
rownames(cs_loc_rain_na) <- cs_loc98_na$REP_ID
cs_loc_rain <- rbind(na.omit(cs_loc_rain),
                     cs_loc_rain_na)
colnames(cs_loc_rain98) <- paste("X",c(1:12), "1998", sep = "_")

cs_loc_rain98_long <- cs_loc_rain98 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  pivot_longer(starts_with("X"), names_to = c("Month", "Year"),
               names_prefix = "X_", names_sep = "_",
               values_to = "rainfall") %>%
  mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "[A-Z]"),"[", 1)),
         rainfall = ifelse(rainfall < 9e20, rainfall, NA)) %>%
  left_join(select(sample_date, SERIES_NUM, DATE = MID_DATE98)) %>%
  mutate(DATE = ifelse(!is.na(DATE), DATE, 
                       round(mean(DATE, na.rm = TRUE))+1),
         Month = as.numeric(Month)) %>%
  mutate(field_season = ifelse(Month <= DATE & Month >= DATE - 3, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_rainfall = mean(rainfall),
            sum_rainfall = sum(rainfall))

# 1990
rain90 <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_199001-199012.nc"
rain <- raster::brick(rain90)
rain[rain > 9e20] <- NA

cs_loc90 <- plot_locations %>%
  filter(YEAR == "y90") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain90 <- raster::extract(rain, cs_loc90)
rownames(cs_loc_rain90) <- cs_loc90$REP_ID
# for reps that don't fall in the raster get average of all values within 2000m
cs_loc90_na <- cs_loc90 %>%
  filter(REP_ID %in% rownames(cs_loc_rain)[is.na(cs_loc_rain[,1])])
cs_loc_rain_na <- raster::extract(rain, cs_loc90_na,
                                  fun = mean, buffer = 2000)
rownames(cs_loc_rain_na) <- cs_loc90_na$REP_ID
cs_loc_rain <- rbind(na.omit(cs_loc_rain),
                     cs_loc_rain_na)
colnames(cs_loc_rain90) <- paste("X",c(1:12), "1990", sep = "_")

cs_loc_rain90_long <- cs_loc_rain90 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  pivot_longer(starts_with("X"), names_to = c("Month", "Year"),
               names_prefix = "X_", names_sep = "_",
               values_to = "rainfall") %>%
  mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "[A-Z]"),"[", 1)),
         rainfall = ifelse(rainfall < 9e20, rainfall, NA)) %>%
  left_join(select(sample_date, SERIES_NUM, DATE = MID_DATE90)) %>%
  mutate(DATE = ifelse(!is.na(DATE), DATE, 
                       round(mean(DATE, na.rm = TRUE))+1),
         Month = as.numeric(Month)) %>%
  mutate(field_season = ifelse(Month <= DATE & Month >= DATE - 3, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_rainfall = mean(rainfall),
            sum_rainfall = sum(rainfall))


# 1978
rain78 <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_197801-197812.nc"
rain <- raster::brick(rain78)
rain[rain > 9e20] <- NA

cs_loc78 <- plot_locations %>%
  filter(YEAR == "y78") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain78 <- raster::extract(rain, cs_loc78)
rownames(cs_loc_rain78) <- cs_loc78$REP_ID
# for reps that don't fall in the raster get average of all values within 2000m
cs_loc78_na <- cs_loc78 %>%
  filter(REP_ID %in% rownames(cs_loc_rain)[is.na(cs_loc_rain[,1])])
cs_loc_rain_na <- raster::extract(rain, cs_loc78_na,
                                  fun = mean, buffer = 2000)
rownames(cs_loc_rain_na) <- cs_loc78_na$REP_ID
cs_loc_rain <- rbind(na.omit(cs_loc_rain),
                     cs_loc_rain_na)
colnames(cs_loc_rain78) <- paste("X",c(1:12), "1978", sep = "_")

cs_loc_rain78_long <- cs_loc_rain78 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  pivot_longer(starts_with("X"), names_to = c("Month", "Year"),
               names_prefix = "X_", names_sep = "_",
               values_to = "rainfall") %>%
  mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "[A-Z]"),"[", 1)),
         rainfall = ifelse(rainfall < 9e20, rainfall, NA)) %>%
  left_join(select(sample_date, SERIES_NUM, DATE = MID_DATE78)) %>%
  mutate(DATE = ifelse(!is.na(DATE), DATE, 
                       round(mean(DATE, na.rm = TRUE))+1),
         Month = as.numeric(Month)) %>%
  mutate(field_season = ifelse(Month <= DATE & Month >= DATE - 3, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_rainfall = mean(rainfall),
            sum_rainfall = sum(rainfall))

#combine years
cs_survey_rainfall <- do.call(rbind, list(cs_loc_rain07_long,
                                          cs_loc_rain98_long,
                                          cs_loc_rain90_long,
                                          cs_loc_rain78_long))
ggplot(cs_survey_rainfall, aes(x = mean_rainfall)) +
  geom_histogram() +
  facet_wrap(~Year) +
  labs(x = "Average monthly rainfall for 4 months pre-survey")

p1 <- ggplot(cs_survey_rainfall, aes(x = mean_rainfall)) +
  geom_histogram() +
  facet_wrap(~Year, ncol = 1, scales = "free_y") +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Average monthly rainfall for 4 months pre-survey")
p2 <- ggplot(cs_survey_rainfall, aes(x = sum_rainfall)) +
  geom_histogram() +
  facet_wrap(~Year, ncol = 1, scales = "free_y") +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Total rainfall for 4 months pre-survey")
p1 + p2

ggplot(cs_survey_rainfall, aes(x = mean_rainfall, y = sum_rainfall)) + 
  geom_point() +
  geom_abline(slope=4, intercept = 0) +
  facet_wrap(~Year)
# basically the same

p1
ggsave("Average monthly rainfall for 4 months pre-survey.png",
       path = "Outputs/Graphs/",width = 12, height = 15, units = "cm")


# Climatic averages
str(plot_locations)
CS19_locs <- plot_locations %>%
  filter(REP_ID %in% VEGETATION_PLOTS_20161819$REP_ID) %>%
  mutate(YEAR = "y19") %>%
  group_by(REP_ID, YEAR, SERIES_NUM) %>%
  summarise_if(is.numeric, mean) %>%
  mutate(ID = NA, REPEAT_PLOT_ID = NA, OS_8_FIG_10M = NA,
         OS_6_FIG_100M = NA, OS_4_FIG_1KM = NA,OS_2_FIG_10KM = NA,) %>%
  select(all_of(colnames(plot_locations)))
allplot_loc <- plot_locations %>%
  rbind(CS19_locs) %>%
  mutate(plot_x = ifelse(!is.na(E_6_FIG_100M), E_6_FIG_100M, E_4_FIG_1KM),
         plot_y = ifelse(!is.na(N_6_FIG_100M), N_6_FIG_100M, N_4_FIG_1KM)) %>%
  select(REP_ID, YEAR, plot_x, plot_y) %>%
  na.omit() %>%
  st_as_sf(coords = c("plot_x","plot_y"), crs = 27700)

clim_rain <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_ann-30y_198101-201012.nc"
rain <- raster::brick(clim_rain)
cs_loc_rain30y <- raster::extract(rain, allplot_loc)
colnames(cs_loc_rain30y) <- "RAIN_8110"
cs_loc_rain30y <- cbind(as.data.frame(allplot_loc),cs_loc_rain30y)
cs_loc_rain30y <- cs_loc_rain30y %>%
  select(-geometry) %>%
  mutate(RAIN_8110 = ifelse(RAIN_8110 < 9e20, RAIN_8110, NA),
         Year = recode(YEAR,
                       "y19" = 2019,
                       "y07" = 2007,
                       "y9899" = 1998,
                       "y90" = 1990,
                       "y78" = 1978))

sum_rain <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_seas-30y_198101-201012.nc"
rain <- raster::brick(sum_rain)
cs_loc_sumrain30y <- raster::extract(rain, allplot_loc)
colnames(cs_loc_sumrain30y) <- c("WIN","SPR","SUM","AUT")
cs_loc_sumrain30y <- cbind(as.data.frame(allplot_loc),cs_loc_sumrain30y)
cs_loc_sumrain30y <- cs_loc_sumrain30y %>%
  select(-geometry) %>%
  mutate(across(WIN:AUT, function(x) ifelse(x < 9e20,x, NA)),
         Year = recode(YEAR,
                       "y19" = 2019,
                       "y07" = 2007,
                       "y9899" = 1998,
                       "y90" = 1990,
                       "y78" = 1978))

cs_rainfall_stats <- full_join(cs_loc_rain30y, cs_loc_sumrain30y) %>%
  full_join(mutate(ungroup(cs_survey_rainfall), Year = as.numeric(Year))) %>%
  mutate(rain_diff = mean_rainfall - SUM/3) %>%
  select(REP_ID, Year, AVER_RAIN_8110 = RAIN_8110, rain_diff)

cs_rainfall_diff <- cs_survey_rainfall %>% ungroup() %>%
  na.omit() %>% filter(Year != 1990) %>%
  filter(grepl("X", REP_ID)) %>%
  select(Year:mean_rainfall) %>%
  pivot_wider(id_cols = REP_ID, names_from = Year, 
              values_from = mean_rainfall,
              names_prefix = "rain") %>%
  mutate(diff7898 = rain1998 - rain1978,
         diff9807 = rain2007 - rain1998) %>%
  select(REP_ID, contains("diff")) %>%
  pivot_longer(contains("diff"), values_to = "fieldseason_rain",
               names_to = "Time_period", names_prefix = "diff") %>%
  mutate(Year = recode(Time_period,
                       "7898" = 1998,
                       "9807" = 2007)) %>%
  filter(grepl("X", REP_ID))

cs_rainfall_averages <-
  full_join(select(cs_loc_rain30y, REP_ID, Year, AVER_RAIN = RAIN_8110),
            select(cs_loc_sumrain30y, REP_ID, Year, AVER_SUM_RAIN = SUM))
  

# CHESS soil moisture ####
nc <- ncdf4::nc_open("~/Shapefiles/CHESS/CHESSLandMonHydEn2007.nc")
nc
soil_moist <- ncdf4::ncvar_get(nc, attributes(nc$var)$names[25])
easting <- ncdf4::ncvar_get(nc, "eastings")
northing <- ncdf4::ncvar_get(nc, "northings")
# take top layer of soil only
soil_moist <- soil_moist[,1,] %>% cbind(easting, northing) 
colnames(soil_moist) <- c(month.abb, "easting","northing")
moist_crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs"

soil_moist <- as.data.frame(soil_moist)
moist <- st_as_sf(soil_moist, coords = c("easting","northing"), crs = 27700)

cs_loc_moist <- st_nearest_feature(cs_loc07, moist)
cs_loc_moist <- as.data.frame(soil_moist)[cs_loc_moist,1:12]
rownames(cs_loc_moist) <- cs_loc07$REP_ID
colnames(cs_loc_moist) <- paste0("X_",1:12)

cs_loc_moist07_long <- cs_loc_moist %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  pivot_longer(starts_with("X_"), names_to = "Month",
               values_to = "moisture", names_prefix = "X_") %>%
  mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "[A-Z]"),"[", 1)),
         Year = 2007, Month = as.numeric(Month)) %>%
  left_join(select(sample_date, SERIES_NUM, DATE = MID__DATE07)) %>%
  mutate(field_season = ifelse(Month <= DATE & Month >= DATE - 3, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_moisture = mean(moisture), month = max(Month))

cs_loc_moist07_sample_month <- cs_loc_moist %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  pivot_longer(starts_with("X_"), names_to = "Month",
               values_to = "moisture", names_prefix = "X_") %>%
  mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "[A-Z]"),"[", 1)),
         Year = 2007, Month = as.numeric(Month)) %>%
  left_join(select(sample_date, SERIES_NUM, DATE = MID__DATE07)) %>%
  filter(Month == DATE) %>%
  group_by(Year, REP_ID) %>%
  summarise(eo_moisture = mean(moisture))




# Atmospheric deposition data ####
# Cumulative N deposition calculation
str(Ndep_avg)

Ndep_avg_cumsum <- Ndep_avg %>%
  mutate(Navg_cumdep78 = rowSums(select(.,gridavg_1970:gridavg_1978)),
         Navg_cumdep90 = rowSums(select(.,gridavg_1970:gridavg_1990)),
         Navg_cumdep98 = rowSums(select(.,gridavg_1970:gridavg_1998)),
         Navg_cumdep07 = rowSums(select(.,gridavg_1970:gridavg_2007)),
         Navg_cumdep18 = rowSums(select(.,gridavg_1970:gridavg_2018)),
         Navg_cumdep07_30 = rowSums(select(.,gridavg_1977:gridavg_2007)),
         Navg_cumdep18_30 = rowSums(select(.,gridavg_1988:gridavg_2018))) %>%
  select(x,y,contains("cumdep"))

Ndep_for_cumsum <- Ndep_for %>%
  mutate(Nfor_cumdep78 = rowSums(select(.,forest_1970:forest_1978)),
         Nfor_cumdep90 = rowSums(select(.,forest_1970:forest_1990)),
         Nfor_cumdep98 = rowSums(select(.,forest_1970:forest_1998)),
         Nfor_cumdep07 = rowSums(select(.,forest_1970:forest_2007)),
         Nfor_cumdep18 = rowSums(select(.,forest_1970:forest_2018)),
         Nfor_cumdep07_30 = rowSums(select(.,forest_1977:forest_2007)),
         Nfor_cumdep18_30 = rowSums(select(.,forest_1988:forest_2018))) %>%
  select(x,y,contains("cumdep"))

Ndep_moo_cumsum <- Ndep_moo %>%
  mutate(Nmoo_cumdep78 = rowSums(select(.,moor_1970:moor_1978)),
         Nmoo_cumdep90 = rowSums(select(.,moor_1970:moor_1990)),
         Nmoo_cumdep98 = rowSums(select(.,moor_1970:moor_1998)),
         Nmoo_cumdep07 = rowSums(select(.,moor_1970:moor_2007)),
         Nmoo_cumdep18 = rowSums(select(.,moor_1970:moor_2018)),
         Nmoo_cumdep07_30 = rowSums(select(.,moor_1977:moor_2007)),
         Nmoo_cumdep18_30 = rowSums(select(.,moor_1988:moor_2018))) %>%
  select(x,y,contains("cumdep"))

Ndep_cumsum <- full_join(Ndep_avg_cumsum, 
                         Ndep_for_cumsum) %>%
  full_join(Ndep_moo_cumsum)

Ndep_cumsum70 <- Ndep_cumsum %>%
  select(-ends_with("_30")) %>%
  melt(id.vars = c("x","y"), variable.factor = FALSE, 
       value.name = "Ndep") %>%
  mutate(Year = recode(substring(variable,12,13),
                       "78" = 1978,
                       "90" = 1990,
                       "98" = 1998,
                       "07" = 2007,
                       "18" = 2018),
         Habitat = recode(substring(variable, 2,4),
                          "avg" = "gridavg",
                          "for" = "forest",
                          "moo" = "moor")) %>% select(-variable)

# Rate of change for Sdep
str(Sdep_avg)
Sdep_avg_change <- Sdep_avg %>%
  mutate(Savg_change78 = gridavg_1978 - gridavg_1970,
         Savg_change90 = gridavg_1990 - gridavg_1970,
         Savg_change98 = gridavg_1998 - gridavg_1970,
         Savg_change07 = gridavg_2007 - gridavg_1970,
         Savg_change18 = gridavg_2018 - gridavg_1970,
         Savg_change07_30 = gridavg_2007 - gridavg_1977,
         Savg_change18_30 = gridavg_2018 - gridavg_1988) %>%
  select(x,y,contains("change"))

Sdep_for_change <- Sdep_for %>%
  mutate(Sfor_change78 = forest_1978 - forest_1970,
         Sfor_change90 = forest_1990 - forest_1970,
         Sfor_change98 = forest_1998 - forest_1970,
         Sfor_change07 = forest_2007 - forest_1970,
         Sfor_change18 = forest_2018 - forest_1970,
         Sfor_change07_30 = forest_2007 - forest_1977,
         Sfor_change18_30 = forest_2018 - forest_1988) %>%
  select(x,y,contains("change"))

Sdep_moo_change <- Sdep_moo %>%
  mutate(Smoo_change78 = moor_1978 - moor_1970,
         Smoo_change90 = moor_1990 - moor_1970,
         Smoo_change98 = moor_1998 - moor_1970,
         Smoo_change07 = moor_2007 - moor_1970,
         Smoo_change18 = moor_2018 - moor_1970,
         Smoo_change07_30 = moor_2007 - moor_1977,
         Smoo_change18_30 = moor_2018 - moor_1988) %>%
  select(x,y,contains("change"))

Sdep_change <- full_join(Sdep_avg_change, 
                         Sdep_for_change) %>%
  full_join(Sdep_moo_change)

Sdep_change70 <- Sdep_change %>%
  select(-ends_with("_30")) %>%
  melt(id.vars = c("x","y"), variable.factor = FALSE, 
       value.name = "Sdep") %>%
  mutate(Year = recode(substring(variable,12,13),
                       "78" = 1978,
                       "90" = 1990,
                       "98" = 1998,
                       "07" = 2007,
                       "18" = 2018),
         Habitat = recode(substring(variable, 2,4),
                          "avg" = "gridavg",
                          "for" = "forest",
                          "moo" = "moor")) %>% select(-variable)

AtmosDep_70 <- full_join(Ndep_cumsum70, Sdep_change70)
AtmosDep_70_nona <- na.omit(AtmosDep_70)

# merge with CS plots
# quick check on NAs in plot locations
table(filter(plot_locations, is.na(E_6_FIG_100M))$E_4_FIG_1KM)
table(filter(plot_locations, is.na(N_6_FIG_100M))$N_4_FIG_1KM)
# if there is an NA in the 100m measurement can use the 1km measurement (1km
# measurements has 0s if there are 100m measurements)
plot_locs19 <- plot_locations %>%
  filter(REP_ID %in% VEGETATION_PLOTS_20161819$REP_ID) %>%
  mutate(YEAR = "y19", ID = NA, E_10_FIG_1M = NA, N_10_FIG_1M = NA,
         REPEAT_PLOT_ID = NA, OS_8_FIG_10M = NA, OS_6_FIG_100M = NA) %>%
  group_by(REP_ID) %>%
  mutate(E_6_FIG_100M = mean(E_6_FIG_100M),
         N_6_FIG_100M = mean(N_6_FIG_100M)) %>% ungroup() %>%
  unique() 
janitor::get_dupes(plot_locs19, REP_ID) #0 dupes

# get habitat information for every plot - if no info use gridavg
CS_habs <- BH_comb %>% 
  mutate(Habitat = ifelse(BH %in% c(1,2), "forest", "moor")) %>%
  mutate(Habitat = replace_na(Habitat, "gridavg")) %>% 
  select(REP_ID, Year, Habitat) %>%
  unique()
CS_habs_dupes <- janitor::get_dupes(CS_habs, REP_ID, Year) %>%
  filter(Habitat == "forest") %>% select(-dupe_count)
CS_habs <- CS_habs %>%
  filter(!paste0(REP_ID, Year) %in% 
           paste0(CS_habs_dupes$REP_ID, CS_habs_dupes$Year)) %>%
  rbind(CS_habs_dupes)
janitor::get_dupes(CS_habs, REP_ID, Year)

# get locations of every plot
CS_plot_atdep <- plot_locations %>% 
  rbind(plot_locs19) %>%
  mutate(plot_x = ifelse(!is.na(E_6_FIG_100M), E_6_FIG_100M, E_4_FIG_1KM),
         plot_y = ifelse(!is.na(N_6_FIG_100M), N_6_FIG_100M, N_4_FIG_1KM)) %>%
  select(plot_x, plot_y, REP_ID, Year = YEAR) %>%
  mutate(Year = recode(Year, 
                       "y19" = 2018,
                       "y07" = 2007,
                       "y9899" = 1998,
                       "y90" = 1990,
                       "y78" = 1978)) %>%
  left_join(CS_habs) %>%
  mutate(Habitat = replace_na(Habitat, "gridavg"))

# match to atmospheric deposition data
dep_x <- AtmosDep_70_nona$x
dep_y <- AtmosDep_70_nona$y

for(i in 1:nrow(CS_plot_atdep)) {
  CS_plot_atdep[i,"x"] <- dep_x[which.min(abs(dep_x - CS_plot_atdep$plot_x[i]))]
  CS_plot_atdep[i,"y"] <- dep_y[which.min(abs(dep_y - CS_plot_atdep$plot_y[i]))]
}

CS_plot_atdep <- left_join(CS_plot_atdep, AtmosDep_70_nona)
summary(CS_plot_atdep)
