###------------------------------------------------------------##
## Process data inputs for state-space models
###------------------------------------------------------------##

library(stringi)
library(zoo)
source(here::here('config.R'))

## This script requires that `site` and `pft` be defined. If calling from 
## run_ssm.R, these will be passed in from the shell script.

## Directory of input data files for this site. 
## The site directory must include the following files:
## tree_level_df.csv, plot_level_df.csv, obs_matrix.csv, year_matrix.csv, 
## size_init.csv, and cwd_mean.csv

## Routes us to the correct input data file based on site alias
indir <- file.path('Input_data',indir_names[[site]])

## Full tree list for specified site
allData <- read.table(here::here(indir, "tree_level_df.csv"), header = TRUE, sep = 
                        ",",row.names=1)

## Subset full tree list for trees of specified PFT
treeData <- allData[allData$PFT == pft,]

## If running yellow pine at W. Lake Tahoe (UPFU), ignore the treated unit that 
## doesn't have a control counterpart
if(indir_names[[site]]=='WLakeTahoe' & pft == 'yellow_pine'){
  treeData <- treeData[treeData$PlotID != 'MCK 13-3 T', ]
}

## Plot-level attributes
allPlots <- read.csv(here::here(indir, 'plot_level_df.csv'))

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
obs_all <- read.table(here::here(indir, 'obs_matrix.csv'), header = TRUE, 
                      sep = ',', row.names = 1)
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
year_obs_all <- read.table(here::here(indir, 'year_matrix.csv'), header=TRUE, 
                           sep = ',', row.names = 1)
year_obs <- year_obs_all[row.names(treeData),]

## Read in linearly interpolated DBHs
sizes_interp_all <- read.csv(here::here(indir,'size_init.csv'),row.names = 1)

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
cwd <- read.table(here::here(indir, "cwd_mean.csv"), header=TRUE, 
                  sep = ',', row.names=1)
cwd_pft <- cwd[row.names(treeData),]

## Scale CWD values
cwd_mean = mean(as.matrix(cwd_pft),na.rm = TRUE)
cwd_stdv <- sd(as.matrix(cwd_pft), na.rm = TRUE)
cwd_scaled <- (cwd_pft - cwd_mean)/cwd_stdv