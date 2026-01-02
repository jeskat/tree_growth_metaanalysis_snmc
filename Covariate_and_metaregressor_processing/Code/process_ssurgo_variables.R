### This workflow is adapted from Elissa Chasen (2020), "Cleaning SSURGO data
### for mapping ecological processes," downloaded from 
### https://rpubs.com/emchasen/SSURGOcleaning. 

library(FedData)
library(dplyr)

### Function to remove columns that include only null values
not_all_na <- function(x) any(!is.na(x)) #for data cleaning

### Get data from each survey area 
ca607 <- get_ssurgo(template = 'CA607', label = 'CA607') # LaTour
ca693 <- get_ssurgo(template = 'CA693', label = 'CA693') # W. Lake Tahoe
ca724 <- get_ssurgo(template = 'CA724', label = 'CA724') # Blodgett
ca731 <- get_ssurgo(template = 'CA731', label = 'CA731') # STEF
ca750 <- get_ssurgo(template = 'CA750', label = 'CA750') # Teakettle
ca792 <- get_ssurgo(template = 'CA792', label = 'CA792') # Sequoia and Tharp's Creek

### Combine chorizon, component, and map unit tables from each area
areas_list <- list(ca607, ca693, ca724, ca731, ca750, ca792)
chorizon <- vector(mode = 'list', length = 6)
component <- vector(mode = 'list', length = 6)
mapunit <- vector(mode = 'list', length = 6)

for(a in 1:6){
  area <- areas_list[[a]]
  chorizon[[a]] <- area$tabular$chorizon
  component[[a]] <- area$tabular$component
  mapunit[[a]] <- area$tabular$mapunit
}

chorizon <- do.call(rbind, chorizon)
component <- do.call(rbind, component)
mapunit <- do.call(rbind, mapunit)

### Remove empty columns
chorizon <- chorizon %>%
  select_if(not_all_na)
component <- component %>% 
  select_if(not_all_na)
mapunit <- mapunit %>%
  select_if(not_all_na)

### Determine soil depth (deepeest horizon bottom of each component)
depth <- chorizon %>% 
  group_by(cokey) %>%
  summarise(total_depth = max(hzdepb.r))

### Remove horizons that start below 50 cm
chorizon <- chorizon %>% 
  filter(hzdept.r < 50) %>%
  droplevels()

### Take depth-weighted average across components
chorizon <- chorizon %>%
  mutate(thick = ifelse(hzdepb.r > 50, 50-hzdept.r,
                        hzdepb.r - hzdept.r)) %>%
  group_by(cokey) %>% 
  summarise(ksat = round(weighted.mean(ksat.r, thick, na.rm = TRUE), 2),
            cec = round(weighted.mean(cec7.r, thick, na.rm = TRUE), 2),
            awc = round(weighted.mean(awc.r, thick, na.rm = TRUE), 2))

## Add depth information
chorizon <- left_join(chorizon, depth, by = 'cokey')


### Subset columns of interest from component and mapunit data
component <- component %>% 
  dplyr::select(c(comppct.r, compname, majcompflag, mukey, cokey))

mapunit <- mapunit %>% 
  dplyr::select(c(musym, muname, muacres, mukey))

### Merge soil characteristics with component keys to get component-level data
component_horizon <- left_join(component, chorizon, 
                               by = c('cokey'))

### Merge component-level data with map unit information
full_soil <- left_join(component_horizon, mapunit, by = c('mukey'))

### Replace commas with underscores before exporting to CSV
full_soil <- full_soil %>%
  mutate(muname = gsub(', ', '_', muname))

### Save to disk
write.csv(full_soil, 
          here::here('Covariate_and_metaregressor_processing/Processed_data/Soil/all_unit_soils.csv'), 
          row.names = FALSE)
