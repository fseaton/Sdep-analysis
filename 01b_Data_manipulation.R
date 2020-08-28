# Script for calculating Ellenberg scores for the most recent survey years X plots
library(dplyr)

str(VEGETATION_PLOT_SP_161819)
table(VEGETATION_PLOT_SP_161819$PLOT_TYPE)
unique(VEGETATION_PLOT_SP_161819[VEGETATION_PLOT_SP_161819$PLOT_TYPE=="XX","REP_ID"])
table(VEGETATION_PLOT_SP_161819[VEGETATION_PLOT_SP_161819$PLOT_TYPE=="X","PLOTYEAR"])

str(SPECIES_LIB_TRAITS)
filter(SPECIES_LIB_CODES, COLUMN_NAME == "GROWTH_FORM")

# AVC data manipulation
hab07 <- select(CS07_IBD, REP_ID = REP_ID07, AVC07) %>%
  unique()
hab90 <- select(CS90_IBD, REP_ID = REP_ID90, AVC90) %>%
  unique()
hab98 <- select(CS98_IBD, REP_ID = REP_ID98, AVC98) %>%
  unique()
hab78 <- select(CS78_IBD, REP_ID = REP_ID78, AVC78) %>%
  unique()

# create combined AVC variable, if 07 has AVC use that otherwise use 98 then 78.
# There are only 3 sites with no AVC data and I can't see how to get theirs as
# they don't appear in 2016/19.
hab <- full_join(hab07, hab98) %>%
  full_join(hab90) %>%
  full_join(hab78) %>%
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


## Ellenberg scores ####
CS18_ELL <- filter(VEGETATION_PLOT_SP_161819, PLOT_TYPE %in% c("X","XX")) %>%
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
  left_join(select(VEGETATION_PLOTS_2007, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
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
                   REP_ID, Year, BRC_NUMBER, TOTAL_COVER = COVER)) %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   starts_with("EBER"))) %>%
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% # set 0 values to NA
  group_by(Year, REP_ID) %>%
  summarise(across(starts_with("EBER"), mean, na.rm = TRUE)) %>%
  rename_with(~gsub("EBERG","SM_",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA))

X_Ell_whole <- CS07_SP %>%
  left_join(select(VEGETATION_PLOTS_2007, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
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

p1 <- ggplot(filter(X_Ell_comp, Year != 1978), aes(x = WH_R, y = SM_R)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc) +
  labs(x = bquote("Ellenberg R 400m"^2~"plot"),
       y = bquote("Ellenberg R 4m"^2~"plot"))

ggplot(filter(X_Ell_comp, Year != 1978), aes(x = WH_N, y = SM_N)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~AVC_desc)

ggplot(filter(X_Ell_comp, Year != 1978), aes(x = R_diff)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~AVC_desc)

ggplot(filter(X_Ell_comp, Year != 1978), aes(x = N_diff)) +
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
  left_join(select(VEGETATION_PLOTS_2007, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
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
  left_join(select(VEGETATION_PLOTS_2007, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
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

library(patchwork)
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

# Quick fix for UK19_PH values that aren't matching to the veg plots
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

