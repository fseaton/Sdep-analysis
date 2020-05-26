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

# veg data
CS78_IBD <- fread("../../../CS/Data/IBD78.csv")
CS98_IBD <- fread("../../../CS/Data/IBD98.csv")
CS07_IBD <- fread("../../../CS/Data/IBD07.csv")

# metadata
BHPH_NAMES <- read.csv("../../../CS/Data/BHPH_names.csv")
