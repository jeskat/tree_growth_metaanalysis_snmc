###------------------------------------------------------------##
## Process data inputs for state-space models
###------------------------------------------------------------##

library(stringi)
library(zoo)
library(dplyr)
library(tidyr)
source(here::here('State_space_growth_models/Code/ssm_config.R'))

## This script requires that `site` and `pft` be defined. If calling from 
## run_ssm.R, these will be passed in from the shell script.

## Directory of input data files for this site. 
## The site directory must include the following files:
## tree_attrs.csv, unit_attrs.csv, dbh_tree_obs.csv, year_tree_obs.csv, 
## pft_df.csv and cwd_mean.csv

## Routes us to the correct input data file based on site alias
indir <- file.path('State_space_growth_models/','Input_data',indir_names[[site]])

## Full tree list for specified site
allData <- read.table(here::here(indir, "tree_attrs.csv"), header = TRUE, sep = 
                        ",",row.names='TreeID')

## Subset full tree list for trees of specified PFT
treeData <- allData[allData$PFT == pft,]

## If running yellow pine at W. Lake Tahoe (UPFU), ignore the treated unit that 
## doesn't have a control counterpart
if(indir_names[[site]]=='WLakeTahoe' & pft == 'yellow_pine'){
  treeData <- treeData[treeData$PlotID != 'MCK 13-3 T', ]
}

## Unit-level attributes
allPlots <- read.csv(here::here(indir, 'unit_attrs.csv'))

## reorder rows by plotID
allPlots <- allPlots[order(allPlots$PlotID),]
rownames(allPlots) <- (allPlots$PlotID)

## Subset plots that have >=1 tree of the specified PFT
plotData <- subset(allPlots, rownames(allPlots) %in% unique(treeData$PlotID))

## Make plotID a factor
treeData$fPlot <- factor(treeData$PlotID)

## Make treatment a factor
plotData$fTreatment <- factor(plotData$treatment)

## Make None the reference level for treatment
plotData$fTreatment <- relevel(plotData$fTreatment, 'None')


## Load DBH observations for all trees and then subset to the desired pft
obs_all <- read.table(here::here(indir, 'dbh_tree_obs.csv'), header = TRUE, 
                      sep = ',', row.names = 'TreeID')
obs_pft <- obs_all[row.names(treeData),]

## Get the constants needed for scaling and centering transformed DBH values
obs_mean = mean(as.matrix(obs_pft),na.rm = TRUE)
obs_stdv =  sd(as.matrix(obs_pft), na.rm = TRUE)

obs_sq <- obs_pft**2
obs_sq_mean <- mean(as.matrix(obs_sq), na.rm = TRUE)
obs_sq_sd <- sd(as.matrix(obs_sq), na.rm = TRUE)

obs_ln <- log(obs_pft)
obs_ln_mean <- mean(as.matrix(obs_ln), na.rm = TRUE)
obs_ln_stdv <- sd(as.matrix(obs_ln), na.rm = TRUE)


## Load years corresponding to each dbh observation, and subset to desired pft
year_obs_all <- read.table(here::here(indir, 'year_tree_obs.csv'), header=TRUE, 
                           sep = ',', row.names = 'TreeID')
year_obs <- year_obs_all[row.names(treeData),]

### Get initial DBH values for each tree and each year
## If the initial DBH values have already been created, use them.
if(file.exists(here::here(indir, 'initDBH_tree_year.csv'))){
  sizes_interp_all <- read.table(here::here(indir, 'initDBH_tree_year.csv'), header=TRUE, 
                                 sep = ',', row.names = 'TreeID')
}else{
  ## Otherwise, create initial DBH values from the processed inventory data.
  ## Read full dataframe with inventory years
  pft_df <- read.table(here::here(indir, "pft_df.csv"), header = TRUE, sep = 
                         ",")
  
  ## Determine the full range of years and the first year of the study
  # The first_yr should be the minimum 'Year' across the entire study.
  first_yr <- min(pft_df$Year)
  
  # 2. Create the full time series for every tree and calculate the Timeseries_Year
  size_init <- pft_df %>%
    # For each tree (TreeID)...
    group_by(TreeID) %>%
    # ...generate a complete sequence of years from min year to max year for that tree
    tidyr::complete(Year = seq(min(Year), max(Year), by = 1)) %>%
    ungroup() %>%
    # Calculate the Timeseries_Year column (1, 2, 3, ...)
    mutate(
      Timeseries_Year = Year - first_yr + 1
    ) %>%
    # Select and arrange the columns
    select(TreeID, Year, Timeseries_Year, DBH)
  
  # 3. Linearly Interpolate Missing DBH Values
  # The linear interpolation in R is handled by the na.approx function from the 'zoo' package
  size_init <- size_init %>%
    group_by(TreeID) %>%
    # Interpolate the missing DBH values (NA) within each tree's group
    # na.approx() performs linear interpolation
    mutate(
      DBH = na.approx(DBH, na.rm = FALSE)
    ) %>%
    ungroup()
  
  # 4. Pivot the Data Table to  wide format
  # The output is a data frame where each row is a tree and columns are Timeseries_Year
  sizes_interp_all <- size_init %>%
    # Use pivot_wider to convert from long to wide format
    # The names_from column becomes the new column headers
    # The values_from column provides the values for the new cells
    pivot_wider(
      id_cols = TreeID,
      names_from = Timeseries_Year,
      values_from = DBH
    ) %>%
    # Ensure the Timeseries_Year columns are ordered numerically, though pivot_wider often handles this.
    arrange(TreeID) %>%
    tibble::column_to_rownames(var = "TreeID")
}

## Subset to the desired PFT
sizes_interp_pft <- sizes_interp_all[row.names(treeData),]

## Create initial values for all latent sizes (DBHs) in all years
nTree <- nrow(obs_pft)
nYr <- max(allData$end_year) - min(allData$start_year) +1

size_init_pft <- matrix(nrow = nTree, 
                        ncol = nYr)
for(i in 1:nTree){
  ## The first latent DBH is selected from a normal distribution centered on the observed DBH
  size_init_pft[i, treeData$start_year[i]] <- abs(rnorm(1, mean = obs_pft[i, 1], 
                                                        sd = 0.1)) 
  for(yr in (treeData$start_year[i]+1):treeData$end_year[i]){
    ## For subsequent years, the initial value is the max of the previous latent 
    ## DBH + 0.001 or the absolute random value from a normal distribution 
    ## centered on the observed diameter. This ensures latend dDBH is positive.
    size_init_pft[i, yr] = max(size_init_pft[i, yr-1] + 0.001,
                               abs(rnorm(1, mean = sizes_interp_pft[i, yr], 
                                         sd = 0.1)))
  }
}

## Make the ddbh matrices that correspond to the size matrices
matching_ddbh_init <- matrix(nrow = nTree, ncol = nYr)
matching_ddbh_init[, 1] <- NA
matching_ddbh_init[,2:nYr] <- size_init_pft[,2:nYr] - size_init_pft[,1:(nYr-1)]

## Make a set of initial sizes to fill for only the first column
matching_size_init <- matrix(nrow = nTree, ncol = nYr)

## First column of size is stochastic; the rest of size is deterministic
for(i in 1:nTree){
  matching_size_init[i, treeData$start_year[i]] <- size_init_pft[i, treeData$start_year[i]]
}

## Load CWD data and subset for PFT
cwd <- read.table(here::here(indir, "cwd_tree_year.csv"), header=TRUE, 
                  sep = ',', row.names='TreeID')
cwd_pft <- cwd[row.names(treeData),]

## Scale CWD values
cwd_mean = mean(as.matrix(cwd_pft),na.rm = TRUE)
cwd_stdv <- sd(as.matrix(cwd_pft), na.rm = TRUE)
cwd_scaled <- (cwd_pft - cwd_mean)/cwd_stdv