# Code for reading in data 
library(data.table)
library(readxl)

# Laptop
# Atmospheric data
Sdep_avg <- fread("../CBED_1970-2018_SN_dep/totalCBED_Sdep_gridavg_1970-2018_5km_kgS_ha.csv")
Sdep_for <- fread("../CBED_1970-2018_SN_dep/totalCBED_Sdep_forest_1970-2018_5km_kgS_ha.csv")
Sdep_moo <- fread("../CBED_1970-2018_SN_dep/totalCBED_Sdep_moor_1970-2018_5km_kgS_ha.csv")

Ndep_avg <- fread("../CBED_1970-2018_SN_dep/totalCBED_Ndep_gridavg_1970-2018_5km_kgN_ha.csv")
Ndep_for <- fread("../CBED_1970-2018_SN_dep/totalCBED_Ndep_forest_1970-2018_5km_kgN_ha.csv")
Ndep_moo <- fread("../CBED_1970-2018_SN_dep/totalCBED_Ndep_moor_1970-2018_5km_kgN_ha.csv")

# Soils data
CS78_PH <- fread("../../../CS/Data/CS1978_PH_LOI_DATA.csv")
CS98_PH <- fread("../../../CS/Data/CS1998_PH_LOI_DATA.csv")
CS07_PH <- fread("../../../CS/Data/CS2007_PH_LOI_DATA.csv")
CS16_PH <- fread("../../../CS/Data/CS2016_PH_LOI_DATA.csv")
UK19_PH <- fread("../../UK19_PHLOI.csv")

CS_tier4 <- fread("../../../CS/Data/CS_SOILS_TIER_4_DATA.csv")
UK19_WET <- read_excel("../../Copy of UKSCAPE Field Survey 2019 CHEMICAL DARv3.xls",
                       sheet = "DATA", skip = 3, na = "NA")

# veg data
CS78_IBD <- fread("../../../CS/Data/IBD78.csv")
CS90_IBD <- fread("../../../CS/Data/IBD90.csv")
CS98_IBD <- fread("../../../CS/Data/IBD98.csv")
CS07_IBD <- fread("../../../CS/Data/IBD07.csv")
GM16_IBD <- fread("../../../CS/Data/GMEP_VEG_IBD_XNEST1.csv")

# Plant species records
CS19_SP <- fread("../../../CS/Data/VEGETATION_PLOT_SP_161819.csv")
CS07_SP <- fread("../../../CS/Data/VEGETATION_PLOT_SPECIES_2007.csv")
CS98_SP <- fread("../../../CS/Data/CS9899_SPECIESDATA.csv")
CS90_SP <- fread("../../../CS/Data/CS90_SPECIESDATA.csv")
CS78_SP <- fread("../../../CS/Data/CS78_SPECIESDATA.csv")

# metadata
BHPH_NAMES <- read.csv("../../../CS/Data/BHPH_names.csv")
CS07_PLOTS <- read.csv("../../../CS/Data/VEGETATION_PLOTS_2007.csv")
landclass_dat <- read.csv("../../../CS/Data/CS_landclass_dat.csv")
SPECIES_LIB_TRAITS <- fread("../../../CS/Data/SPECIES_LIB_TRAITS.csv")
SPECIES_LIB_CODES <- fread("../../../CS/Data/SPECIES_LIB_CODES.csv")

# Desktop
library(RODBC)
## get data on species characteristics from CS database
pwds <- read.csv("Outputs/pwd.csv")
channel3 <- odbcConnect("FEGEN", uid="csgeo", pwd=pwds[pwds$uid=="csgeo","pwd"], believeNRows=FALSE)
sqlTables(channel3, schema="CSVEG")
channelg <- odbcConnect("tbb", uid="AXISII_WP3_VIEWER", 
                        pwd=pwds[pwds$uid=="AXISII_WP3_VIEWER","pwd"], believeNRows=FALSE)


# derived IBD files
CS78_IBD <- sqlFetch(channel3, "CSVEG.IBD78")
CS90_IBD <- sqlFetch(channel3, "CSVEG.IBD90")
CS98_IBD <- sqlFetch(channel3, "CSVEG.IBD98")
CS07_IBD <- sqlFetch(channel3, "CSVEG.IBD07")
GM16_IBD <- sqlFetch(channelg, "GMEP_DERIVED.WP6_IBD_YRS1234")

# metadata
SPECIES_LIB_TRAITS <- sqlFetch(channel3, "CSVEG.LUS_SP_LIB_AND_TRAITS")
SPECIES_LIB_CODES <- sqlFetch(channel3, "CSVEG.LUS_SP_LIB_CODES_NEW")

# recent plant data
CS19_SP <- sqlFetch(channel3, "CSVEG.VEGETATION_PLOT_SP_161819") 
CS78_SP <- sqlFetch(channel3, "CSVEG.CS78_SPECIESDATA")
CS90_SP <- sqlFetch(channel3, "CSVEG.CS90_SPECIESDATA")
CS98_SP <- sqlFetch(channel3, "CSVEG.CS9899_SPECIESDATA")
CS07_SP <- sqlFetch(channel3, "CSVEG.VEGETATION_PLOT_SPECIES_2007")
VEGETATION_PLOTS_20161819 <- sqlFetch(channel3, "CSVEG.VEGETATION_PLOTS_20161819")
VEGETATION_PLOTS_2007 <- sqlFetch(channel3, "CSVEG.VEGETATION_PLOTS_2007")

# soils data
channel2 <- odbcConnect("MWA", uid="masq", pwd=pwds[pwds$uid=="masq","pwd"], believeNRows=FALSE)
sqlTables(channel2, schema="DB_MASQ")

# pH data
CS16_PH <- sqlFetch(channel2, "DB_MASQ.CS2016_PH_LOI_DATA")
CS07_PH <- sqlFetch(channel2, "DB_MASQ.CS2007_PH_LOI_DATA")
CS78_PH <- sqlFetch(channel2, "DB_MASQ.CS1978_PH_LOI_DATA")
CS98_PH <- sqlFetch(channel2, "DB_MASQ.CS1998_PH_LOI_DATA")
UK19_PH <- read.csv("Outputs/UK19_PHLOI.csv")

# tier 4 data
CS_tier4 <- sqlFetch(channel2, "DB_MASQ.CS_SOILS_TIER_4_DATA")

# landclass data
channel <- odbcConnect("MWA", uid = "lus_user", pwd = pwds[pwds$uid=="lus_user","pwd"], 
                       case = "nochange", believeNRows = FALSE)
sqlTables(channel, schema = "DB_CSSQUARES")
landclass_dat <- sqlFetch(channel, "DB_CSSQUARES.SQUARES_FILE_ALL_LC")
