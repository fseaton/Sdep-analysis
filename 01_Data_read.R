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

# Plant species recordS
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
# atmospheric deposition data
Sdep_avg <- fread("N:/UKScape/CBED_1970-2018_SN_dep/totalCBED_Sdep_gridavg_1970-2018_5km_kgS_ha.csv")
Sdep_for <- fread("N:/UKScape/CBED_1970-2018_SN_dep/totalCBED_Sdep_forest_1970-2018_5km_kgS_ha.csv")
Sdep_moo <- fread("N:/UKScape/CBED_1970-2018_SN_dep/totalCBED_Sdep_moor_1970-2018_5km_kgS_ha.csv")

Ndep_avg <- fread("N:/UKScape/CBED_1970-2018_SN_dep/totalCBED_Ndep_gridavg_1970-2018_5km_kgN_ha.csv")
Ndep_for <- fread("N:/UKScape/CBED_1970-2018_SN_dep/totalCBED_Ndep_forest_1970-2018_5km_kgN_ha.csv")
Ndep_moo <- fread("N:/UKScape/CBED_1970-2018_SN_dep/totalCBED_Ndep_moor_1970-2018_5km_kgN_ha.csv")

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

# soils data
MWA_masq <- dbConnect(odbc(), "MWA", UID="masq", 
                      PWD=pwds[pwds$uid=="masq","pwd"])

# pH data
CS16_PH <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS2016_PH_LOI_DATA"))
CS07_PH <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS2007_PH_LOI_DATA"))
CS78_PH <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS1978_PH_LOI_DATA"))
CS98_PH <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS1998_PH_LOI_DATA"))
UK19_PH <- read.csv("Outputs/UK19_PHLOI.csv")

CS07_CN <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS2007_CN_DATA"))
CS98_CN <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS1998_CN_DATA"))

# tier 4 data
CS_tier4 <- dbReadTable(MWA_masq, SQL("DB_MASQ.CS_SOILS_TIER_4_DATA"))

# landclass data
MWA_lus <- dbConnect(odbc(), "MWA", UID="lus_user", 
                           PWD=pwds[pwds$uid=="lus_user","pwd"])

landclass_dat <- dbReadTable(MWA_lus, SQL("DB_CSSQUARES.SQUARES_FILE_ALL_LC"))
plot_locations <- dbReadTable(MWA_lus, SQL("DB_CSSQUARES.SQUARE_AND_PLOT_LOCATIONS"))
plot_dates <- dbReadTable(MWA_lus, SQL("DB_CSSQUARES.DATEVARS_78_84_90_98_07"))

# CEH laptop data
UK19_WET <- read_excel("C:/Users/fseaton/Documents/CS/Copy of UKSCAPE Field Survey 2019 CHEMICAL DARv3.xls",
                       sheet = "DATA", skip = 3, na = "NA")
