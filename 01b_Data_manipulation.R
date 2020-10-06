# Script for various data manipulations
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(sf)
library(data.table)

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

Soil_missingreps78 <- CS78_PH$REP_ID[!CS78_PH$REP_ID %in% CS_REPEAT_PLOTS$Y78]
Soil_missingreps98 <- CS98_PH$REP_ID[!CS98_PH$REP_ID %in% CS_REPEAT_PLOTS$Y9899]
Soil_missingreps07 <- CS07_PH$REP_ID[!CS07_PH$REP_ID %in% CS_REPEAT_PLOTS$Y07]

# now correct repeat plots for 2019 data and missing soils plots
CS_REP_ID <- CS_REPEAT_PLOTS %>% #filter(CS_REPEAT_PLOTS, AMALG_PTYPE == "X") %>%
  mutate(Y07 = ifelse(!is.na(Y07), Y07,
                      ifelse(Y9899 %in% Soil_missingreps07, Y9899, NA)),
         Y78 = ifelse(!is.na(Y78), Y78, 
                      ifelse(Y90 %in% Soil_missingreps78, Y90, NA)),
         Y9899 = ifelse(!is.na(Y9899), Y9899,
                        ifelse(Y90 %in% Soil_missingreps98, Y90, NA))) %>%
  mutate(Y19 = Y07)

CS_REP_ID_LONG <- CS_REP_ID %>%
  select(REPEAT_PLOT_ID, Y78:Y19) %>%
  pivot_longer(Y78:Y19, names_to = "Year",
               values_to = "REP_ID") %>%
  na.omit() %>%
  mutate(Year = recode(Year,
                       "Y78" = 1978,
                       "Y90" = 1990,
                       "Y9899" = 1998,
                       "Y07" = 2007,
                       "Y19" = 2019))

rm(list = ls(pattern = "^Soil_missingreps"))

# AVC data manipulation ####
hab07 <- select(CS07_IBD, REP_ID = REP_ID07, AVC07) %>%
  unique() %>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, Y07), 
            by = c("REP_ID" = "Y07")) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID)) %>%
  select(-REPEAT_PLOT_ID)
hab98 <- select(CS98_IBD, REP_ID = REP_ID98, AVC98) %>%
  unique()%>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, Y9899), 
            by = c("REP_ID" = "Y9899")) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID))%>%
  select(-REPEAT_PLOT_ID)
hab90 <- select(CS90_IBD, REP_ID = REP_ID90, AVC90) %>%
  unique()%>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, Y90), 
            by = c("REP_ID" = "Y90")) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID))%>%
  select(-REPEAT_PLOT_ID)
hab78 <- select(CS78_IBD, REP_ID = REP_ID78, AVC78) %>%
  unique()%>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, Y78), 
            by = c("REP_ID" = "Y78")) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID))%>%
  select(-REPEAT_PLOT_ID)

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
# Ellenberg for inner 2x2m square
X_Ell_inner <- CS07_SP %>%
  left_join(select(CS07_PLOTS, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
  filter(PLOT_TYPE == "X" & NEST_LEVEL < 2)  %>%
  mutate(Year = 2007) %>%
  select(REP_ID, Year, BRC_NUMBER, FIRST_COVER) %>%
  full_join(select(mutate(filter(CS19_SP, 
                                 PLOT_TYPE %in% c("X","XX") & NEST_LEVEL < 2), Year = 2019),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER)) %>%
  full_join(select(mutate(filter(CS98_SP, grepl("X", REP_ID) & NEST_LEVEL < 2), Year = 1998),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER)) %>%
  full_join(select(mutate(filter(CS90_SP, grepl("X", REP_ID) & QUADRAT_NEST < 2), Year = 1990),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER = C1)) %>%
  full_join(select(mutate(filter(CS78_SP, grepl("X", REP_ID) & QUADRAT_NEST < 2), Year = 1978),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER = COVER)) %>%
  na.omit() %>%
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
  na.omit() %>%
  left_join(unique(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                          starts_with("EBER")))) %>%
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% # set 0 values to NA
  group_by(Year, REP_ID) %>%
  summarise(across(starts_with("EBER"), mean, na.rm = TRUE)) %>%
  rename_with(~gsub("EBERG","WH_",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA))

# weighted mean Ellenberg 
X_wEll_inner <- CS07_SP %>%
  left_join(select(CS07_PLOTS, VEG_PLOTS_ID, REP_ID, PLOT_TYPE)) %>%
  filter(PLOT_TYPE == "X" & NEST_LEVEL < 2)  %>%
  mutate(Year = 2007) %>%
  select(REP_ID, Year, BRC_NUMBER, FIRST_COVER) %>%
  full_join(select(mutate(filter(CS19_SP, 
                                 PLOT_TYPE %in% c("X","XX") & NEST_LEVEL < 2), Year = 2019),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER)) %>%
  full_join(select(mutate(filter(CS98_SP, grepl("X", REP_ID) & NEST_LEVEL < 2), Year = 1998),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER)) %>%
  full_join(select(mutate(filter(CS90_SP, grepl("X", REP_ID) & QUADRAT_NEST < 2), Year = 1990),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER = C1)) %>%
  full_join(select(mutate(filter(CS78_SP, grepl("X", REP_ID) & QUADRAT_NEST < 2), Year = 1978),
                   REP_ID, Year, BRC_NUMBER, FIRST_COVER = COVER)) %>%
  na.omit() %>%
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
  na.omit() %>%
  left_join(unique(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                          starts_with("EBER")))) %>%
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% 
  group_by(Year, REP_ID) %>%
  summarise(across(starts_with("EBER"), weighted.mean, w = TOTAL_COVER, na.rm = TRUE)) %>%
  rename_with(~gsub("EBERG","WH_",.x)) %>%
  mutate_all(function(x) ifelse(!is.nan(x), x, NA))

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
  ) %>%
  ungroup() %>%
  mutate(REP_ID = gsub("XX","X",REP_ID)) %>%
  left_join(CS_REP_ID_LONG) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID)) %>%
  select(-REPEAT_PLOT_ID)



# pH data ####
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
  full_join(select(UK19_PH, REP_ID, PH_2019 = PH_DIW, PHC_2019 = PH_CACL)) %>%
  unique() %>%
  mutate(PH_2019 = ifelse(!is.na(PH_2019), PH_2019, PH_2016),
         PHC_2019 = ifelse(!is.na(PHC_2019), PHC_2019, PHC_2016)) %>%
  select(-ends_with("2016"))
str(PH)
summary(PH)
mice::md.pattern(PH)
janitor::get_dupes(PH, REP_ID)

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
         diff7819 = PH_2019 - PH_1978,
         diff9807 = PH_2007 - PH_1998,
         diff9819 = PH_2019 - PH_1998,
         diff0719 = PH_2019 - PH_2007) 
summary(PH)

PH_diff_long <- PH %>%
  select(REP_ID, starts_with("diff")) %>%
  pivot_longer(starts_with("diff"),
               values_to = "pH",
               values_drop_na = TRUE) %>%
  mutate(name = as.factor(name)) %>%
  mutate(name = forcats::fct_inorder(name))

# select only most recent change and convert into wide format for plotting
PH_Diff_wide <- select(PH, REP_ID, diff0719) %>%
  na.omit() %>%
  left_join(select(CS07_PLOTS, REP_ID, POINT_X, POINT_Y))



PH_QA_diff <- data.frame(
  Time_period = c("7898","9807","0719"),
  PH_DIW_SE = c(0.357,0.303,0.187),
  PH_CACL2_SE = c(NA,NA,0.065),
  rain_diff_sd = 22
)


# LOI ####
LOI <- full_join(select(CS78_PH, REP_ID, LOI_1978 = LOI1978),
                 select(CS98_PH, REP_ID, LOI_1998 = LOI2000)) %>%
  full_join(select(CS07_PH, REP_ID, LOI_2007 = LOI2007)) %>%
  full_join(select(CS16_PH, REP_ID, LOI_2016 = LOI)) %>%
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
                         names_transform = list(Year = as.integer)) %>%
  mutate(Year = ifelse(Year == 2016, 2019, Year))
LOI_long


# Soil moisture ####
UK19_WET <- UK19_WET %>%
  mutate(REP_ID = paste0(`Square number`,`X plot`)) 

MOISTURE <- CS_tier4 %>%
  mutate(REP_ID = ifelse(!is.na(REP_ID07), REP_ID07,
                         ifelse(!is.na(REP_ID98), REP_ID98,
                                ifelse(!is.na(REP_ID78), REP_ID78, NA)))) %>%
  select(REP_ID, MOISTURE_CONTENT_07, MOISTURE_CONTENT_98) %>%
  full_join(select(UK19_WET, REP_ID, MOISTURE_CONTENT_19 = `g_water/wet_weight_of_soil`)) %>%
  # mutate(VWC_19 = 100*VWC_19) %>%
  pivot_longer(starts_with("MOISTURE"), names_to = "variable", 
               values_to = "Moisture", values_drop_na = TRUE) %>%
  mutate(Year = ifelse(variable == "MOISTURE_CONTENT_07", 2007,
                       ifelse(variable == "MOISTURE_CONTENT_98", 1998,
                              ifelse(variable == "MOISTURE_CONTENT_19", 2019, NA)))) %>%
  select(-variable) %>% na.omit() %>% unique() %>%
  left_join(CS_REP_ID_LONG) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID)) %>%
  select(-REPEAT_PLOT_ID)


# Spatial data manipulation ####
# rainfall - monthly totals ####
str(plot_dates)
sample_date <- plot_dates %>%
  select(SERIES_NUM, starts_with("MID_")) %>%
  full_join(
    select(VEGETATION_PLOTS_20161819, SERIES_NUM = SQUARE, 
           MID_DATE19 = COMPLETED_DATE) %>%
      unique()) %>%
  mutate(across(starts_with("MID_"), lubridate::ymd)) %>%
  mutate(across(starts_with("MID_"), lubridate::month))

# 2018
rain18 <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_201801-201812.nc"
rain <- raster::brick(rain18)
rain[rain > 9e20] <- NA

cs_loc18 <- VEGETATION_PLOTS_20161819 %>%
  filter(PLOTYEAR == 2018) %>%
  select(REP_ID, POINT_X, POINT_Y) %>%
  # have 3 duplicated plots so taking average of the points given
  # 377D4, 431Y1, 912R1
  group_by(REP_ID) %>%
  summarise(POINT_X = mean(POINT_X), POINT_Y = mean(POINT_Y)) %>%
  na.omit() %>%
  st_as_sf(coords = c("POINT_X","POINT_Y"), crs = 27700)

cs_loc_rain18 <- raster::extract(rain, cs_loc18)
rownames(cs_loc_rain18) <- cs_loc18$REP_ID
# all reps fall in raster
colnames(cs_loc_rain18) <- paste("X",c(1:12), "2018", sep = "_")


cs_loc_rain18_long <- cs_loc_rain18 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  pivot_longer(starts_with("X"), names_to = c("Month", "Year"),
               names_prefix = "X_", names_sep = "_",
               values_to = "rainfall") %>%
  mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "[A-Z]"),"[", 1)),
         rainfall = ifelse(rainfall < 9e20, rainfall, NA)) %>%
  left_join(select(sample_date, SERIES_NUM, DATE = MID_DATE19)) %>%
  mutate(DATE = ifelse(!is.na(DATE), DATE, 
                       round(mean(DATE, na.rm = TRUE))+1),
         Month = as.numeric(Month)) %>%
  mutate(field_season = ifelse(Month <= DATE & Month > DATE - 4, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_rainfall = mean(rainfall),
            sum_rainfall = sum(rainfall)) 


# 2016
rain16 <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_201601-201612.nc"
rain <- raster::brick(rain16)
rain[rain > 9e20] <- NA

cs_loc16 <- VEGETATION_PLOTS_20161819 %>%
  filter(PLOTYEAR == 2016) %>%
  select(REP_ID, POINT_X, POINT_Y) %>%
  # have 3 duplicated plots so taking average of the points given
  # 377D4, 431Y1, 912R1
  group_by(REP_ID) %>%
  summarise(POINT_X = mean(POINT_X), POINT_Y = mean(POINT_Y)) %>%
  na.omit() %>%
  st_as_sf(coords = c("POINT_X","POINT_Y"), crs = 27700)

cs_loc_rain16 <- raster::extract(rain, cs_loc16)
rownames(cs_loc_rain16) <- cs_loc16$REP_ID
# for reps that don't fall in the raster get average of all values within 2000m
cs_loc16_na <- cs_loc16 %>%
  filter(REP_ID %in% rownames(cs_loc_rain16)[is.na(cs_loc_rain16[,1])])
cs_loc_rain16_na <- raster::extract(rain, cs_loc16_na,
                                    fun = mean, buffer = 2000)
rownames(cs_loc_rain16_na) <- cs_loc16_na$REP_ID
cs_loc_rain16 <- rbind(na.omit(cs_loc_rain16),
                       cs_loc_rain16_na)
colnames(cs_loc_rain16) <- paste("X",c(1:12), "2016", sep = "_")


cs_loc_rain16_long <- cs_loc_rain16 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  pivot_longer(starts_with("X"), names_to = c("Month", "Year"),
               names_prefix = "X_", names_sep = "_",
               values_to = "rainfall") %>%
  mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "[A-Z]"),"[", 1)),
         rainfall = ifelse(rainfall < 9e20, rainfall, NA)) %>%
  left_join(select(sample_date, SERIES_NUM, DATE = MID_DATE19)) %>%
  mutate(DATE = ifelse(!is.na(DATE), DATE, 
                       round(mean(DATE, na.rm = TRUE))+1),
         Month = as.numeric(Month)) %>%
  mutate(field_season = ifelse(Month <= DATE & Month > DATE - 4, 1,0)) %>%
  filter(field_season == 1) %>%
  group_by(Year, REP_ID) %>%
  summarise(mean_rainfall = mean(rainfall),
            sum_rainfall = sum(rainfall)) 


# 2007
rain07 <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_mon_200701-200712.nc"
rain <- raster::brick(rain07)
rain[rain > 9e20] <- NA

cs_loc07 <- plot_locations %>%
  filter(YEAR == "y07") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain07 <- raster::extract(rain, cs_loc07)
rownames(cs_loc_rain07) <- cs_loc07$REP_ID
# for reps that don't fall in the raster get average of all values within 2000m
cs_loc07_na <- cs_loc07 %>%
  filter(REP_ID %in% rownames(cs_loc_rain07)[is.na(cs_loc_rain07[,1])])
cs_loc_rain07_na <- raster::extract(rain, cs_loc07_na,
                                  fun = mean, buffer = 2000)
rownames(cs_loc_rain07_na) <- cs_loc07_na$REP_ID
cs_loc_rain07 <- rbind(na.omit(cs_loc_rain07),
                     cs_loc_rain07_na)
colnames(cs_loc_rain07) <- paste("X",c(1:12), "2007", sep = "_")


cs_loc_rain07_long <- cs_loc_rain07 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  pivot_longer(starts_with("X"), names_to = c("Month", "Year"),
               names_prefix = "X_", names_sep = "_",
               values_to = "rainfall") %>%
  mutate(SERIES_NUM = as.numeric(sapply(strsplit(REP_ID, "[A-Z]"),"[", 1)),
         rainfall = ifelse(rainfall < 9e20, rainfall, NA)) %>%
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

cs_loc98 <- plot_locations %>%
  filter(YEAR == "y9899") %>%
  dplyr::select(REP_ID, E_10_FIG_1M, N_10_FIG_1M) %>%
  na.omit() %>%
  st_as_sf(coords = c("E_10_FIG_1M","N_10_FIG_1M"), crs = 27700)

cs_loc_rain98 <- raster::extract(rain, cs_loc98)
rownames(cs_loc_rain98) <- cs_loc98$REP_ID
# for reps that don't fall in the raster get average of all values within 2000m
cs_loc98_na <- cs_loc98 %>%
  filter(REP_ID %in% rownames(cs_loc_rain98)[is.na(cs_loc_rain98[,1])])
cs_loc_rain_na <- raster::extract(rain, cs_loc98_na,
                                  fun = mean, buffer = 2000)
rownames(cs_loc_rain_na) <- cs_loc98_na$REP_ID
cs_loc_rain98 <- rbind(na.omit(cs_loc_rain98),
                       cs_loc_rain_na)
colnames(cs_loc_rain98) <- paste("X",c(1:12), "1998", sep = "_")

cs_loc_rain98_long <- cs_loc_rain98 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  # mutate(REP_ID = substring(REP_ID, 2)) %>%
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
  filter(REP_ID %in% rownames(cs_loc_rain90)[is.na(cs_loc_rain90[,1])])
cs_loc_rain_na <- raster::extract(rain, cs_loc90_na,
                                  fun = mean, buffer = 2000)
rownames(cs_loc_rain_na) <- cs_loc90_na$REP_ID
cs_loc_rain90 <- rbind(na.omit(cs_loc_rain90),
                       cs_loc_rain_na)
colnames(cs_loc_rain90) <- paste("X",c(1:12), "1990", sep = "_")

cs_loc_rain90_long <- cs_loc_rain90 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  # mutate(REP_ID = substring(REP_ID, 2)) %>%
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
  filter(REP_ID %in% rownames(cs_loc_rain78)[is.na(cs_loc_rain78[,1])])
cs_loc_rain_na <- raster::extract(rain, cs_loc78_na,
                                  fun = mean, buffer = 2000)
rownames(cs_loc_rain_na) <- cs_loc78_na$REP_ID
cs_loc_rain78 <- rbind(na.omit(cs_loc_rain78),
                       cs_loc_rain_na)
colnames(cs_loc_rain78) <- paste("X",c(1:12), "1978", sep = "_")

cs_loc_rain78_long <- cs_loc_rain78 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("REP_ID") %>%
  # mutate(REP_ID = substring(REP_ID, 2)) %>%
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
cs_survey_rainfall <- do.call(rbind, list(cs_loc_rain18_long,
                                          cs_loc_rain16_long,
                                          cs_loc_rain07_long,
                                          cs_loc_rain98_long,
                                          cs_loc_rain90_long,
                                          cs_loc_rain78_long))
cs_survey_rainfall <- CS_REP_ID_LONG %>%
  mutate(Year = as.character(Year)) %>%
  na.omit() %>%
  right_join(mutate(cs_survey_rainfall, 
                    Year = ifelse(Year %in% c("2016","2018"), "2019",Year))) %>%
  mutate(REP_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID)) %>%
  select(-REPEAT_PLOT_ID)


# Climatic averages
str(plot_locations)
CS19_locs <- VEGETATION_PLOTS_20161819 %>%
  select(SERIES_NUM = SQUARE, REP_ID, YEAR = PLOTYEAR,
         E_10_FIG_1M = POINT_X, N_10_FIG_1M = POINT_Y)%>%
  left_join(select(CS_REP_ID, REPEAT_PLOT_ID, REP_ID = Y19)) %>%
  mutate(REPEAT_PLOT_ID = ifelse(!is.na(REPEAT_PLOT_ID), REPEAT_PLOT_ID, REP_ID),
         ID = NA, REPEAT_PLOT_ID = NA, OS_8_FIG_10M = NA,
         OS_6_FIG_100M = NA, OS_4_FIG_1KM = NA,OS_2_FIG_10KM = NA,
         E_6_FIG_100M = round(E_10_FIG_1M, -2), 
         N_6_FIG_100M = round(N_10_FIG_1M, -2), 
         E_4_FIG_1KM = round(E_10_FIG_1M, -3), 
         N_4_FIG_1KM = round(N_10_FIG_1M, -3), 
         E_2_FIG_10KM = round(E_10_FIG_1M, -4), 
         N_2_FIG_10KM = round(N_10_FIG_1M, -4)) %>%
  select(all_of(colnames(plot_locations)))
  
allplot_loc <- plot_locations %>%
  filter(REP_ID != "SQ_BL") %>%
  rbind(CS19_locs) %>%
  mutate(plot_x = ifelse(!is.na(E_6_FIG_100M), E_6_FIG_100M, E_4_FIG_1KM),
         plot_y = ifelse(!is.na(N_6_FIG_100M), N_6_FIG_100M, N_4_FIG_1KM)) %>%
  select(REPEAT_PLOT_ID, plot_x, plot_y) %>%
  na.omit() %>% unique() %>%
  group_by(REPEAT_PLOT_ID) %>%
  summarise(plot_x = mean(plot_x), plot_y = mean(plot_y)) %>%
  st_as_sf(coords = c("plot_x","plot_y"), crs = 27700)

clim_rain <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_ann-30y_198101-201012.nc"
rain <- raster::brick(clim_rain)
rain[rain>9e20] <- NA
cs_loc_rain30y <- raster::extract(rain, allplot_loc)
rownames(cs_loc_rain30y) <- allplot_loc$REPEAT_PLOT_ID
# for reps that don't fall in the raster get average of all values within 2000m
cs_loc_rain30y_na <- allplot_loc %>%
  filter(REPEAT_PLOT_ID %in% rownames(cs_loc_rain30y)[is.na(cs_loc_rain30y[,1])])
cs_loc_rain_na <- matrix(raster::extract(rain, cs_loc_rain30y_na,
                                         fun = mean, buffer = 2000))
rownames(cs_loc_rain_na) <- cs_loc_rain30y_na$REPEAT_PLOT_ID
cs_loc_rain30y <- rbind(na.omit(cs_loc_rain30y),
                        cs_loc_rain_na)
colnames(cs_loc_rain30y) <- "RAIN_8110"
cs_loc_rain30y <- cbind(as.data.frame(allplot_loc),cs_loc_rain30y)
cs_loc_rain30y <- cs_loc_rain30y %>%
  select(-geometry) 

sum_rain <- "~/Shapefiles/HadUK-Grid/rainfall_hadukgrid_uk_1km_seas-30y_198101-201012.nc"
rain <- raster::brick(sum_rain)
rain[rain>9e20] <- NA
cs_loc_sumrain30y <- raster::extract(rain, allplot_loc)
rownames(cs_loc_sumrain30y) <- allplot_loc$REPEAT_PLOT_ID
# for reps that don't fall in the raster get average of all values within 2000m
cs_loc_sumrain30y_na <- allplot_loc %>%
  filter(REPEAT_PLOT_ID %in% rownames(cs_loc_sumrain30y)[is.na(cs_loc_sumrain30y[,1])])
cs_loc_rain_na <- raster::extract(rain, cs_loc_sumrain30y_na,
                                  fun = mean, buffer = 2000)
rownames(cs_loc_rain_na) <- cs_loc_sumrain30y_na$REPEAT_PLOT_ID
cs_loc_sumrain30y <- rbind(na.omit(cs_loc_sumrain30y),
                           cs_loc_rain_na)
colnames(cs_loc_sumrain30y) <- c("WIN","SPR","SUM","AUT")
cs_loc_sumrain30y <- cbind(as.data.frame(allplot_loc),cs_loc_sumrain30y)
cs_loc_sumrain30y <- cs_loc_sumrain30y %>%
  select(-geometry) 

cs_rainfall_stats <- full_join(cs_loc_rain30y, cs_loc_sumrain30y) %>%
  rename(REP_ID = REPEAT_PLOT_ID) %>%
  full_join(mutate(ungroup(cs_survey_rainfall), Year = as.numeric(Year))) %>%
  mutate(rain_diff = mean_rainfall - SUM/3) %>%
  select(REP_ID, Year, AVER_RAIN_8110 = RAIN_8110, rain_diff) 

cs_rainfall_diff <- cs_survey_rainfall %>% ungroup() %>%
  na.omit() %>% filter(Year != 1990) %>%
  select(Year:mean_rainfall) %>%
  pivot_wider(id_cols = REP_ID, names_from = Year, 
              values_from = mean_rainfall,
              names_prefix = "rain") %>%
  mutate(diff7898 = rain1998 - rain1978,
         diff9807 = rain2007 - rain1998,
         diff0719 = rain2019 - rain2007) %>%
  select(REP_ID, contains("diff")) %>%
  pivot_longer(contains("diff"), values_to = "fieldseason_rain",
               names_to = "Time_period", names_prefix = "diff",
               values_drop_na = TRUE) %>%
  mutate(Year = recode(Time_period,
                       "7898" = 1998,
                       "9807" = 2007,
                       "0719" = 2019)) 

cs_rainfall_averages <-
  full_join(select(cs_loc_rain30y, REP_ID = REPEAT_PLOT_ID, AVER_RAIN = RAIN_8110),
            select(cs_loc_sumrain30y, REP_ID = REPEAT_PLOT_ID, AVER_SUM_RAIN = SUM)) 



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
# get habitat information for every plot - if no info use gridavg
CS_habs <- BH_comb_nodupes %>% 
  mutate(Habitat = ifelse(BH %in% c(1,2), "forest", "moor")) %>%
  mutate(Habitat = replace_na(Habitat, "gridavg")) %>% 
  select(REP_ID, Year, Habitat) %>%
  unique()
janitor::get_dupes(CS_habs, REP_ID, Year)

# get locations of every plot
CS_plot_atdep <- data.frame(REPEAT_PLOT_ID = allplot_loc$REPEAT_PLOT_ID) %>%
  cbind(st_coordinates(allplot_loc)) %>%
  left_join(CS_REP_ID_LONG) %>%
  select(plot_x = X, plot_y = Y, REP_ID = REPEAT_PLOT_ID, Year) %>%
  mutate(Year = ifelse(Year == 2019, 2018, Year)) %>%
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

rm(list=ls(pattern = "Sdep|Ndep|AtmosDep"))
rm(list=ls(pattern="cs_loc"))
rm(list=ls(pattern="rain[0-9]{2}$"))
rm(rain)