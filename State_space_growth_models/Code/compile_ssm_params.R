## This script loops through all SSM models and compiles the final parameter 
## estimates into one dataframe.  

source(here::here('config.R'))
source(here::here('Code/clean_functions.R'))

## Define file paths for list of complete models
complete_models <- readRDS(here::here('Outputs', model_list_fn))

# Instantiate empty list to hold parameters from each model
alldata <- vector('list', length = length(complete_models))

for(i in 1:length(complete_models)){
  ## Get list of sub-directories associated with model
  ## Change paths to "runs_to_concatenate" lists as needed. 
  if(complete_models[[i]][['chain']] == 'c1'){
    runs <- readRDS(here::here('Outputs', runs_c1_fn))
  }else if(complete_models[[i]][['chain']] == 'c2'){
    runs <- readRDS(here::here('Outputs', runs_c2_fn))
  }
  
  s <- complete_models[[i]][['site']]
  m <- complete_models[[i]][['model']]
  p <- complete_models[[i]][['pft']]
  
  outdirs <- runs[[m]][[s]][[p]]
  
  ## Read the parameter summary from the final segement of the run
  df <- read.csv(file.path(mcmc_results_dir, outdirs[length(outdirs)], 
                           paste0('params_', s, '_', p, '_', m, '.csv')))
  
  ## Append to list
  alldata[[i]] <- addSitePFTModCols(df, s, p, m)
  }

## List to dataframe  
new <- do.call(rbind, alldata)

## Add variance
new$variance <- new$sd**2

## Save csv
fn <- paste0('allParams_', Sys.Date(), '.csv')
write.csv(new, here::here('Outputs', 'Parameter_summaries', fn), row.names = FALSE)

  