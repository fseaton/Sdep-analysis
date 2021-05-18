# Code for reading in data 
library(data.table)
library(readxl)
library(dplyr)

# atmospheric deposition data
Sdep_year <- fread("~/Data/NS_Dep_CS_Data/CS_tot_S_dep_1970_2017_kgS_ha.csv")
Ndep_year <- fread("~/Data/NS_Dep_CS_Data/CS_tot_N_dep_1970_2017_kgN_ha.csv")

# atmospheric deposition data at 5km
Sdep_avg <- fread("~/Data/CBED_1970-2018_SN_dep/totalCBED_Sdep_gridavg_1970-2018_5km_kgS_ha.csv")
Sdep_for <- fread("~/Data/CBED_1970-2018_SN_dep/totalCBED_Sdep_forest_1970-2018_5km_kgS_ha.csv")
Sdep_moo <- fread("~/Data/CBED_1970-2018_SN_dep/totalCBED_Sdep_moor_1970-2018_5km_kgS_ha.csv")

Ndep_avg <- fread("~/Data/CBED_1970-2018_SN_dep/totalCBED_Ndep_gridavg_1970-2018_5km_kgN_ha.csv")
Ndep_for <- fread("~/Data/CBED_1970-2018_SN_dep/totalCBED_Ndep_forest_1970-2018_5km_kgN_ha.csv")
Ndep_moo <- fread("~/Data/CBED_1970-2018_SN_dep/totalCBED_Ndep_moor_1970-2018_5km_kgN_ha.csv")

## get data on species characteristics from CS database
pwds <- read.csv("Outputs/pwd.csv")
library(odbc)
library(DBI)
FEGEN <- dbConnect(odbc(), "FEGEN", UID="csgeo", PWD=pwds[pwds$uid=="csgeo","pwd"])

dbListTables(FEGEN, schema = "CSVEG")

channelg <- dbConnect(odbc(), "tbb", UID="AXISII_WP3_VIEWER", 
                        PWD=pwds[pwds$uid=="AXISII_WP3_VIEWER","pwd"])

# derived IBD files
CS78_IBD <- dbReadTable(FEGEN, SQL("CSVEG.IBD78"))
CS90_IBD <- dbReadTable(FEGEN, SQL("CSVEG.IBD90"))
CS98_IBD <- dbReadTable(FEGEN, SQL("CSVEG.IBD98"))
CS07_IBD <- dbReadTable(FEGEN, SQL("CSVEG.IBD07"))
GM16_IBD <- dbReadTable(channelg, SQL("GMEP_DERIVED.WP6_IBD_YRS1234"))

# metadata
SPECIES_LIB_TRAITS <- dbReadTable(FEGEN, SQL("CSVEG.LUS_SP_LIB_AND_TRAITS"))
SPECIES_LIB_CODES <- dbReadTable(FEGEN, SQL("CSVEG.LUS_SP_LIB_CODES_NEW"))
BHPH_NAMES <- dbReadTable(FEGEN, SQL("CSVEG.BHPH_NAMES"))

# plant species records
CS19_SP <- dbReadTable(FEGEN, SQL("CSVEG.VEGETATION_PLOT_SP_161819"))
CS78_SP <- dbReadTable(FEGEN, SQL("CSVEG.CS78_SPECIESDATA"))
CS90_SP <- dbReadTable(FEGEN, SQL("CSVEG.CS90_SPECIESDATA"))
CS98_SP <- dbReadTable(FEGEN, SQL("CSVEG.CS9899_SPECIESDATA"))
CS07_SP <- dbReadTable(FEGEN, SQL("CSVEG.VEGETATION_PLOT_SPECIES_2007"))

# plot records
VEGETATION_PLOTS_20161819 <- dbReadTable(FEGEN, SQL("CSVEG.VEGETATION_PLOTS_20161819"))
CS07_PLOTS <- dbReadTable(FEGEN, SQL("CSVEG.VEGETATION_PLOTS_2007"))
CS98_PLOTS <- dbReadTable(FEGEN, SQL("CSVEG.CS9899_QUADS_DESCRIPTION"))
CS90_PLOTS <- dbReadTable(FEGEN, SQL("CSVEG.CS90_QUADS_DESCRIPTION"))
CS78_PLOTS <- dbReadTable(FEGEN, SQL("CSVEG.CS78_QUADS_DESCRIPTION"))

# Repeat plot IDs
CS_REPEAT_PLOTS <- read.csv("Outputs/REPEAT_PLOTS_v4.csv")

# plant QA
CS_VEG_QA <- haven::read_sas("Outputs/qa_cs_all_yrs.sas7bdat") %>%
  haven::zap_formats()
CS98_VEG_QA <- haven::read_sas("Outputs/cs9899_speciesdata.sas7bdat") %>%
  haven::zap_formats()
CS16_VEG_QA <- haven::read_sas("Outputs/qa_2016_all_species_data.sas7bdat") %>%
  haven::zap_formats()
CS19_VEG_QA <- haven::read_sas("Outputs/qa_2019_species_data.sas7bdat") %>%
  haven::zap_formats()

# soils data
MWA_masq <- dbConnect(odbc(), "MWA", UID="masq", 
                      PWD=pwds[pwds$uid=="masq","pwd"])

# pH data
CS16_PH <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS2016_PH_LOI_DATA"))
CS07_PH <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS2007_PH_LOI_DATA"))
CS78_PH <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS1978_PH_LOI_DATA"))
CS98_PH <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS1998_PH_LOI_DATA"))
UK19_PH <- read.csv("Outputs/UK19_PHLOI.csv")
UK20_PH <- read.csv("Outputs/UKCEHCS2020v1.csv", na.strings = -9999.9999)

# Total C and N data ####
UK19_CN <- read_excel("Outputs/UKSCAPE Field Survey 2019 CHEMICAL DARv5.xls",
                      sheet = "DATA clean", na = c("","NA",-9999)) %>%
  select(SQUARE = `Square number`, X_PLOT =  `X plot`, 
         C_PERCENT = `C-total (%)...25`, N_PERCENT = `N-total(%)`,
         CN_RATIO = `C:N_ratio`) %>%
  mutate(REP_ID = paste0(SQUARE, X_PLOT)) %>% select(-X_PLOT)
GM16_CN <- dbReadTable(channelg, SQL("GMEP_SOILS.LANC_CHEM"))
CS07_CN <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS2007_CN_DATA"))
CS98_CN <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS1998_CN_DATA"))

CS07_MINN <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS2007_MIN_N_DATA"))

# pH QA data
CS78_PH_QA <- read.csv("Outputs/W1971_remeasurement_pH.csv")
# CS98_PH_QA <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS1998_PH_REANALYSED_DATA"))
CS98_PH_QA <- read_excel("Outputs/Copy of Archive CS2000 soils pH remeasurements.xls",
                         sheet = "For database")
colnames(CS98_PH_QA) <- c("SQUARE_NUM", "PLOT_TYPE", "REP_NUM","PH_DIW_QA",
                          "PH_CACL2_QA")
CS07_PH_QA <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS2007_PH_LOI_QA_DUPLICATES"))
UK19_PH_QA <- read_excel("Outputs/UKSCAPE Field Survey 2019 CHEMICAL DARv5.xls",
                         sheet = "pH and EC", skip = 9, na = c("","NA"))
colnames(UK19_PH_QA) <- c("SQUARE","X_PLOT","PH_DIW","PH_CACL2",
                          "EC","blank","rep_blank","PH_DIW_QA",
                          "PH_CACL2_QA","EC_QA")

# tier 4 data
CS_tier4 <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS_SOILS_TIER_4_DATA"))

# landclass data
MWA_lus <- dbConnect(odbc(), "MWA", UID="lus_user", 
                           PWD=pwds[pwds$uid=="lus_user","pwd"])

landclass_dat <- dbReadTable(MWA_lus, SQL("DB_CSSQUARES.SQUARES_FILE_ALL_LC"))
plot_locations <- dbReadTable(MWA_lus, SQL("DB_CSSQUARES.SQUARE_AND_PLOT_LOCATIONS"))
plot_dates <- dbReadTable(MWA_lus, SQL("DB_CSSQUARES.DATEVARS_78_84_90_98_07"))

# CEH laptop data
UK19_WET <- read_excel("Outputs/UKSCAPE Field Survey 2019 CHEMICAL DARv5.xls",
                       sheet = "DATA clean", na = c("NA",-9999))
