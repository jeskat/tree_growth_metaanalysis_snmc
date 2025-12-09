###------------------------------------------------------------##
## Post-process MCMCs to get growth outcomes at the unit level
## for each model, given a DBH and CWD. 
###------------------------------------------------------------##


library(readr)
library(stringr)
library('coda')

source(here::here('State_space_growth_models/Code/ssm_config.R'))
source(here::here('State_space_growth_models/Code/ssm_functions.R'))
complete_models <- readRDS(here::here('State_space_growth_models/Outputs', model_list_fn))


### Instantiate lists to hold unit-specific growth outcomes and model var-covar matrices
growthResults <- vector('list', length = 0)
covMatrices <- vector('list', length = 0)

for(mod_num in 1:5){ #(length(complete_models))){
  site <- complete_models[[mod_num]][['site']]
  model <- complete_models[[mod_num]][['model']]
  pft <- complete_models[[mod_num]][['pft']]
  burn_in <- complete_models[[mod_num]][['burn_in']]
  chain <- complete_models[[mod_num]][['chain']]
  
  ## Load input data
  source(here::here('State_space_growth_models/Code/process_inputs.R'))
  
  ## Load MCMC
  all_samples <- load_mcmc(
    site = site,
    pft = pft,
    model = model,
    burn_in = burn_in,
    chain = chain,
    mcmc_results_dir = mcmc_results_dir,
    to_report = to_report
  )
  
  ## Get plot random effects
  plot_effects <- load_random_plot_effects(
    site = site,
    pft = pft,
    model = model, 
    burn_in = burn_in,
    chain = chain,
    mcmc_results_dir = mcmc_results_dir
  )
  
  ### Compile MCMCs for growth in each plot
  plots <- levels(treeData$fPlot)
  
  for(DBH in c(seq(20, 100, 10))){
    for(CWD in seq(-3, 3, 1)){
      plot_outs <- get_reparameterized_plot_growth(CWD = CWD,
                                                   DBH = DBH,
                                                   plot_ids = plots,
                                                   plotData = plotData,
                                                   all_samples = all_samples,
                                                   obs_ln_mean = obs_ln_mean,
                                                   obs_ln_stdv = obs_ln_stdv,
                                                   obs_sq_mean = obs_sq_mean,
                                                   obs_sq_sd = obs_sq_sd)
      ## Append to lists of all models
      growthResults[[paste(site, pft, 'DBH', DBH, 'CWD', CWD, sep = '_')]] <- 
        addSitePFTModCols(plot_outs$long_tbl, site, pft, model)
      
      ## Variance-covariance matrix
      covMatrices[[paste(site, pft, 'DBH', DBH, 'CWD', CWD, sep = '_')]] <- 
        plot_outs$cv_mat
    }
  }
  
}


### Combine growth outcomes from all models into one dataframe
to_save <- do.call(rbind, growthResults)

### Save to disk; file names include specified values for CWD and DBH
fn <- paste0('growthOutcomes_CWD', as.character(CWD), 
             '_DBH', as.character(DBH), '_',
             Sys.Date())

write.csv(to_save, 
          here::here('Outputs', 'Growth_outcomes', paste0(fn, '.csv')), 
          row.names=FALSE)

saveRDS(covMatrices, 
        here::here('Outputs', 'Growth_outcomes', paste0(fn, '.RData')))
