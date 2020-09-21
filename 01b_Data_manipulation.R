# Script for various data manipulations
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(sf)

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

# pH data ####
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


CS78_PH$REP_ID <- paste(CS78_PH$SQUARE_NUM,CS78_PH$REP_NUM, sep = "X")
CS98_PH$REP_ID <- paste(CS98_PH$SQUARE_NUM,CS98_PH$REP_NUM, sep = "X")
CS07_PH$REP_ID <- paste(CS07_PH$SQUARE_NUM,CS07_PH$REP_NUM, sep = "X")
CS16_PH$REP_ID <- paste(CS16_PH$SQUARE_NUM,CS16_PH$REP_NUM, sep = "X")

# wide format pH
PH <- full_join(select(CS78_PH, REP_ID, PH_1978 = PH1978),
                select(CS98_PH, REP_ID, PH_1998 = PHF2000)) %>%
  full_join(select(CS07_PH, REP_ID, PH_2007 = PH2007_IN_WATER, PHC_2007 = PH2007_IN_CACL2)) %>%
  full_join(select(CS16_PH, REP_ID, PH_2016 = PH_DIW, PHC_2016 = PH_CACL2)) %>%
  full_join(select(UK19_PH, REP_ID, PH_2019 = PH_DIW, PHC_2019 = PH_CACL))
str(PH)
summary(PH)

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
PH <- PH %>%
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


# Spatial data manipulation ####
# rainfall - monthly totals
rain07 <- "../Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_200701-200712.nc"
rain <- raster::brick(rain07)

cs_loc07 <- plot_locations %>%
  filter(YEAR == "y07") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain <- raster::extract(rain, cs_loc07)
rownames(cs_loc_rain) <- cs_loc07$REP_ID
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
  mutate(field_season = ifelse(Month <= DATE & Month >= DATE - 3, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_rainfall = mean(rainfall),
            sum_rainfall = sum(rainfall))



# 1998
rain98 <- "../Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_199801-199812.nc"
rain <- raster::brick(rain98)

cs_loc98 <- plot_locations %>%
  filter(YEAR == "y9899") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain98 <- raster::extract(rain, cs_loc98)
rownames(cs_loc_rain98) <- cs_loc98$REP_ID
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
  mutate(field_season = ifelse(Month <= DATE & Month >= DATE - 3, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_rainfall = mean(rainfall),
            sum_rainfall = sum(rainfall))

# 1990
rain90 <- "../Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_199001-199012.nc"
rain <- raster::brick(rain90)

cs_loc90 <- plot_locations %>%
  filter(YEAR == "y90") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain90 <- raster::extract(rain, cs_loc90)
rownames(cs_loc_rain90) <- cs_loc90$REP_ID
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
  mutate(field_season = ifelse(Month <= DATE & Month >= DATE - 3, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_rainfall = mean(rainfall),
            sum_rainfall = sum(rainfall))


# 1978
rain78 <- "../Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_197801-197812.nc"
rain <- raster::brick(rain78)

cs_loc78 <- plot_locations %>%
  filter(YEAR == "y78") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain78 <- raster::extract(rain, cs_loc78)
rownames(cs_loc_rain78) <- cs_loc78$REP_ID
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
allplot_loc <- plot_locations %>%
  select(REP_ID, YEAR, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

clim_rain <- "../Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_ann-30y_198101-201012.nc"
rain <- raster::brick(clim_rain)
cs_loc_rain30y <- raster::extract(rain, allplot_loc)
colnames(cs_loc_rain30y) <- "RAIN_8110"
cs_loc_rain30y <- cbind(as.data.frame(allplot_loc),cs_loc_rain30y)
cs_loc_rain30y <- cs_loc_rain30y %>%
  select(-geometry) %>%
  mutate(RAIN_8110 = ifelse(RAIN_8110 < 9e20, RAIN_8110, NA),
         Year = recode(YEAR,
                       "y07" = 2007,
                       "y9899" = 1998,
                       "y90" = 1990,
                       "y78" = 1978))

sum_rain <- "../Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_seas-30y_198101-201012.nc"
rain <- raster::brick(sum_rain)
cs_loc_sumrain30y <- raster::extract(rain, allplot_loc)
colnames(cs_loc_sumrain30y) <- c("WIN","SPR","SUM","AUT")
cs_loc_sumrain30y <- cbind(as.data.frame(allplot_loc),cs_loc_sumrain30y)
cs_loc_sumrain30y <- cs_loc_sumrain30y %>%
  select(-geometry) %>%
  mutate(across(WIN:AUT, function(x) ifelse(x < 9e20,x, NA)),
         Year = recode(YEAR,
                       "y07" = 2007,
                       "y9899" = 1998,
                       "y90" = 1990,
                       "y78" = 1978))


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
  filter(REP_ID %in% CS19_PL$REP_ID) %>%
  mutate(YEAR = "y19", ID = NA, E_10_FIG_1M = NA, N_10_FIG_1M = NA,
         REPEAT_PLOT_ID = NA, OS_8_FIG_10M = NA, OS_6_FIG_100M = NA) %>%
  group_by(REP_ID) %>%
  mutate(E_6_FIG_100M = mean(E_6_FIG_100M),
         N_6_FIG_100M = mean(N_6_FIG_100M)) %>% ungroup() %>%
  unique() 
janitor::get_dupes(plot_locs19, REP_ID) #0 dupes

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
                       "y78" = 1978)) #%>%
  # mutate(x = AtmosDep_70_nona$x[which.min(abs(AtmosDep_70_nona$x - plot_x))],
  #        y = AtmosDep_70_nona$y[which.min(abs(AtmosDep_70_nona$y - plot_y))]) %>%
  # full_join(AtmosDep_70_nona)

dep_x <- AtmosDep_70_nona$x
dep_y <- AtmosDep_70_nona$y

for(i in 1:nrow(CS_plot_atdep)) {
  CS_plot_atdep[i,"x"] <- dep_x[which.min(abs(dep_x - CS_plot_atdep$plot_x[i]))]
  CS_plot_atdep[i,"y"] <- dep_y[which.min(abs(dep_y - CS_plot_atdep$plot_y[i]))]
}

CS_plot_atdep <- full_join(CS_plot_atdep, AtmosDep_70_nona)
summary(CS_plot_atdep)
