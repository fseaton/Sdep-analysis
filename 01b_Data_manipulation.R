# Script for calculating Ellenberg scores for the most recent survey years X plots
library(dplyr)

str(VEGETATION_PLOT_SP_161819)
table(VEGETATION_PLOT_SP_161819$PLOT_TYPE)
unique(VEGETATION_PLOT_SP_161819[VEGETATION_PLOT_SP_161819$PLOT_TYPE=="XX","REP_ID"])
table(VEGETATION_PLOT_SP_161819[VEGETATION_PLOT_SP_161819$PLOT_TYPE=="X","PLOTYEAR"])

str(SPECIES_LIB_TRAITS)
filter(SPECIES_LIB_CODES, COLUMN_NAME == "GROWTH_FORM")

CS18_ELL <- filter(VEGETATION_PLOT_SP_161819, PLOT_TYPE %in% c("X","XX")) %>%
  mutate(REP_ID = paste0(SQUARE,PLOT_TYPE,PLOT_NUMBER)) %>%
  mutate(REP_ID = gsub("XX","X",REP_ID)) %>%
  left_join(select(SPECIES_LIB_TRAITS, BRC_NUMBER, 
                   starts_with("EBER"),
                   GROWTH_FORM)) %>%
  filter(GROWTH_FORM %in% c("f","fe","g","m","s","ss","w")) %>% # filter to vascular plants
  mutate(across(starts_with("EBER"), na_if, y = 0)) %>% # set 0 values to NA
  group_by(REP_ID) %>%
  summarise(across(starts_with("EBER"), mean, na.rm = TRUE))
summary(CS18_ELL)
