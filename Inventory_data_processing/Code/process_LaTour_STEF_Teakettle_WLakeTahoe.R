###---------------------------------------------------------------------------##
## Run this script to download and save pre-processed state-space model inputs for   
## LaTour, STEF, Teakettle, and W. Lake Tahoe.  
## Data files save to the State_space_growth_models/Input_data/ directory.
###---------------------------------------------------------------------------##


library(stringr)
library(purrr)
library(dplyr)
library(data.table)

## Download data tables from EDI
source(here::here('Inventory_data_processing/Code/download_processed_data_fourSites.R'))

## Create a data directory for each site
sites <- unique(dt1$Site)

for(s in sites){
  dir_nm <- paste0('State_space_growth_models/Input_data/', s)
  if (!dir.exists(dir_nm)) {
    dir.create(dir_nm)
  }
}

# Save tree_attrs.csv for each site

setDT(dt1)
dt1[, fwrite(.SD, 
             sprintf("State_space_growth_models/Input_data/%s/tree_attrs.csv",
                     .BY$Site),
             sep = ","), 
    by = Site]


# Save unit_attrs.csv for each site

setDT(dt2)
dt2[, fwrite(.SD, 
             sprintf("State_space_growth_models/Input_data/%s/unit_attrs.csv",
                     .BY$Site),
             sep = ","), 
    by = Site]


# Save dbh_tree_obs.csv for each site

setDT(dt3)
dt3[, fwrite(.SD, 
             sprintf("State_space_growth_models/Input_data/%s/dbh_tree_obs.csv",
                     .BY$Site),
             sep = ","), 
    by = Site]


# Save year_tree_obs.csv for each site

setDT(dt4)
dt4[, fwrite(.SD, 
             sprintf("State_space_growth_models/Input_data/%s/year_tree_obs.csv",
                     .BY$Site),
             sep = ","), 
    by = Site]

# Save initDBH_tree_year.csv for each site

setDT(dt5)
dt5[, fwrite(.SD, 
             sprintf("State_space_growth_models/Input_data/%s/initDBH_tree_year.csv",
                     .BY$Site),
             sep = ","), 
    by = Site]

# Save cwd_tree_year.csv for each site

setDT(dt6)
dt6[, fwrite(.SD, 
             sprintf("State_space_growth_models/Input_data/%s/cwd_tree_year.csv",
                     .BY$Site),
             sep = ","), 
    by = Site]
