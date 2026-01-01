### ___________________________________________________________________________
### Takes raw climate variables from Terraclimate; shifts calendar year to 
### climate year and aggregates each variable from monthly to annual.
### ___________________________________________________________________________

library(lubridate)
library(terra)
library(sf)

## Load terraclimate data from netcdf files downloaded from the source
cwd <- rast(
  here::here('Covariate_and_metaregressor_processing/Input_data/TerraClimate_datafiles/agg_terraclimate_def_1958_CurrentYear_GLOBE.nc'))
pdsi <- rast(
  here::here('Covariate_and_metaregressor_processing/Input_data/TerraClimate_datafiles/agg_terraclimate_PDSI_1958_CurrentYear_GLOBE.nc'))
tmax <- rast(
  here::here('Covariate_and_metaregressor_processing/Input_data/TerraClimate_datafiles/agg_terraclimate_tmax_1958_CurrentYear_GLOBE.nc'))
ppt <- rast(
  here::here('Covariate_and_metaregressor_processing/Input_data/TerraClimate_datafiles/agg_terraclimate_ppt_1958_CurrentYear_GLOBE.nc'))

### Function to shift calendar year to climate year 
### (climate year is October - September)
calendar_to_water_yr <- function(ds){ ## input is a dataset
  dte <- time(ds)
  m <- as.numeric(format(dte, "%m"))
  
  ## October to December are assigned to following climate year
  i <- (m>=10)
  dte[i] <- dte[i] %m+% years(1)
  time(ds) <- dte
  return(ds)
}

## Shift calendar year to climate year for each climate variable
tmax_h2oyr <- calendar_to_water_yr(tmax)
ppt_h2oyr <- calendar_to_water_yr(ppt)
pdsi_h2oyr <- calendar_to_water_yr(pdsi)
cwd_h2oyr <- calendar_to_water_yr(cwd)

### Subset tmax for summer months and ppt for winter months
winter_months <- c(10, 11, 12, 1, 2) # a climate year runs Oct-Sept
summer_months = c(6, 7, 8, 9)

## Remove non-summer months from tmax data
dte <- time(tmax_h2oyr)
m <- as.numeric(format(dte, "%m"))
i <- (m %in% summer_months)
summer_tmax <- tmax_h2oyr[[i]]

## Remove non-winter months from ppt data
dte <- time(ppt_h2oyr)
m <- as.numeric(format(dte, "%m"))
i <- (m %in% winter_months)
winter_ppt_h2oyr <- ppt_h2oyr[[i]]

## Aggregate monthly data to the climate year 
annual_smr_tmax <- tapp(summer_tmax, 'years', mean, na.rm=TRUE) ## Mean summer maximum temperature
annual_wnt_ppt <- tapp(winter_ppt_h2oyr, 'years', sum, na.rm=TRUE) ## Winter total precipitation
annual_cwd <- tapp(cwd_h2oyr, 'years', sum, na.rm=TRUE) ## Annual climate water deficit
annual_pdsi <- tapp(pdsi_h2oyr, 'years', mean, na.rm=TRUE) ## Mean annual PDSI

## Save to disk
writeCDF(annual_smr_tmax, 
         here::here('Covariate_and_metaregressor_processing/Outputs/tmax_by_waterYr.nc'), 
                    overwrite = TRUE)
writeCDF(annual_wnt_ppt, 
         here::here('Outputs/ppt_by_waterYr.nc'), 
                    overwrite = TRUE)
writeCDF(annual_cwd, 
         here::here('Outputs/cwd_by_waterYr.nc'), 
                    overwrite=TRUE)
writeCDF(annual_pdsi, 
         here::here('Outputs/pdsi_by_waterYr.nc'), 
         overwrite = TRUE)
