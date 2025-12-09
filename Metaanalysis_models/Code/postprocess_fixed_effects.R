###------------------------------------------------------------##
## Post-process MCMCs to get growth, CWD, and log(DBH) effects by
## treatment for each model
###------------------------------------------------------------##

library(readr)
library(stringr)
library('coda')
library(dplyr)

source(here::here('State_space_growth_models/Code/ssm_config.R'))
source(here::here('State_space_growth_models/Code/ssm_functions.R'))
complete_models <- readRDS(here::here('State_space_growth_models/Outputs', model_list_fn))


### Instantiate lists to hold per-treatment unit- and per-site growth outcomes
### and site-level CWD, and log(DBH) effects and their respective 
### variance-covariance matrices for all models
unit_growthResults <- vector('list', length = length(complete_models))
unit_growthCovs <- vector('list', length = 0)
site_growthResults <- vector('list', length = length(complete_models))
site_growthCovs <- vector('list', length = 0)
cwdResults <- vector('list', length = length(complete_models))
cwdCovs <- vector('list', length = 0)
sizeResults <- vector('list', length = length(complete_models))
sizeCovs <- vector('list', length = 0)

for(mod_num in 1:(length(complete_models))){
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
  
  ### Reparameterize posteriors to calculate growth for each treatment
  ## Reparameterize for a 30cm tree
  ref_dbh <- 30
  
  # Identify which column name includes the log(DBH) parameter
  lnD_nm <- get_lnD_colnm(all_samples)
  
  ## Reparameterize intercept
  beta0_prime <- all_samples[, 'int_overall'] +
    all_samples[, lnD_nm] * (log(ref_dbh) - obs_ln_mean)/obs_ln_stdv +
    all_samples[, 'd2_size'] * (ref_dbh^2 - obs_sq_mean)/obs_sq_sd
  
  ## Reparameterize treatment effects 
  all_trts_prime <- list()
  for(trt in unique(plotData$treatment)){
    all_trts_prime[[paste0('treat_mean[', trt, ']')]] <- 
      all_samples[, paste0('treat_mean[', trt, ']')] + 
      all_samples[, paste0('treat_size[', trt, ']')] * (log(ref_dbh) - obs_ln_mean)/obs_ln_stdv
  }
  
  ### Compile MCMCs for growth in each plot
  plots <- levels(treeData$fPlot)
  
  ## Instantiate lists to hold means and standard deviations for each plot
  out <- vector('list', length = 0)
  plts <- vector('numeric', length = length(plots))
  means <- vector('numeric', length = length(plots))
  sds <- vector('numeric', length = length(plots))
  treats <- vector('character', length(plots))
  
  for(p in 1:length(plots)){
    ## Get treatment associated with plot
    trt <- plotData[as.character(plotData$PlotID) == plots[p], 'treatment']
    nm <- paste0('treat_mean[', trt, ']')
    treats[p] <- nm
    plts[p] <- p
    
    ## Growth in plot is intercept + treatment effect + plot random effect
    out[[p]] <- beta0_prime + ## overall intercept (reparameterized for DBH=30)
      all_trts_prime[[nm]] + ## treatment effect (reparameterized for DBH = 30)
      plot_effects[, paste0('int_plot[', p, ']')] ## random plot effect
    
    means[p] <- summary(out[[p]])$statistics['Mean']
    sds[p] <- summary(out[[p]])$statistics['SD']
  }
  
  df <- data.frame(out)
  names(df) <- plots
  
  ## Summary table with plot ID, treatment, and mean and SD of growth
  long_tbl <- data.frame(unit =plots, variable = treats, mean=means, sd = sds)
  
  ## Append to lists of all models
  unit_growthResults[[mod_num]] <- addSitePFTModCols(long_tbl, site, pft, model)
  
  ## Get variance-covariance matrix
  unit_growthCovs[[paste0(site, '_', pft)]] <- cov(df)
  
  ### Get growth outcomes at treatment level for each site instead of the unit level

  ## Growth outcome = intercept + treatment effect
  y_by_trt <- lapply(all_trts_prime, function(x) x + beta0_prime)
  y_tbl <- as.mcmc(data.frame(y_by_trt))
  colnames(y_tbl) <- names(y_by_trt)

  ## Get posterior mean growth in each treatment
  y_effects <- summarize_effect_by_trt(names(y_by_trt), y_tbl)

  ## Get variance-covariance matrix
  y_cov <- cov(y_tbl)

  ## Append to lists of all models
  site_growthResults[[mod_num]] <- addSitePFTModCols(y_effects, site, pft, model)
  site_growthCovs[[paste0(site, '_', pft)]] <- y_cov
  
  
  ### Reparameterize posteriors to get CWD effects for each treatment
  
  ## Get non-null CWD effect and interaction variables
  cwdNames <- colnames(all_samples)[grep('treat_cwd', colnames(all_samples))]
  cwdNames <- cwdNames[cwdNames %in% c('treat_cwd[NA]')==FALSE]
  
  ## Get posterior CWD for each treatment
  cwd_by_trt <- all_samples[, 'betaCWD'] + ## Main effect
    all_samples[, cwdNames] ## Interactions
  
  ## Get posterior mean for CWD in each treatment
  cwd_effects <- summarize_effect_by_trt(cwdNames, cwd_by_trt)
  
  ## Get variance-covariance matrix for CWD in each treatment
  cwd_cov <- cov(cwd_by_trt)
  
  ## Append to lists of all models
  cwdResults[[mod_num]] <- addSitePFTModCols(cwd_effects, site, pft, model)
  cwdCovs[[paste0(site, '_', pft)]] <- cwd_cov
  
  ### Reparameterize posteriors to get log(DBH) effects for each treatment
  
  ## Get non-null log(DBH) effect and interaction variables
  sizeNames <- colnames(all_samples)[grep('treat_size', colnames(all_samples))]
  sizeNames <- sizeNames[sizeNames %in% c('treat_size[NA]')==FALSE]
  
  ## Get posterior log(DBH) effect for each treatment
  size_by_trt <- all_samples[, lnD_nm] + ## Main effect
    all_samples[, sizeNames] ## Interactions
  
  ## Get posterior mean for log(DBH) in each treatment
  size_effects <- summarize_effect_by_trt(sizeNames, size_by_trt)
  
  ## Get variance-covariance matrix
  size_cov <- cov(size_by_trt) 
  
  ## Append to lists of all models
  sizeResults[[mod_num]] <- addSitePFTModCols(size_effects, site, pft, model)
  sizeCovs[[paste0(site, '_', pft)]] <- size_cov

}

### Combine growth outcomes, CWD effects, and log(DBH) effects into one dataframe
to_save <- rbind(do.call(rbind, site_growthResults), 
                 do.call(rbind, cwdResults),
                 do.call(rbind, sizeResults))


write.csv(to_save, 
          here::here(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/fixed_effects_by_trt_', 
                            Sys.Date(), '.csv')), 
          row.names=FALSE)

write.csv(do.call(rbind, unit_growthResults), 
          here::here(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/unit_growth_effects_by_trt_', 
                            Sys.Date(), '.csv')), 
          row.names=FALSE)

write.csv(do.call(rbind, site_growthResults), 
          here::here(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/site_growth_effects_by_trt_', 
                            Sys.Date(), '.csv')), 
          row.names=FALSE)

saveRDS(unit_growthCovs, 
        here::here(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/cov_unit_growth_by_trt_', 
                          Sys.Date(), '.RData')))

saveRDS(site_growthCovs, 
        here::here(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/cov_site_growth_by_trt_', Sys.Date(), '.RData')))

saveRDS(cwdCovs, 
        here::here(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/cov_cwd_by_trt_', Sys.Date(), '.RData')))

saveRDS(sizeCovs, 
        here::here(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/cov_size_by_trt_', Sys.Date(), '.RData')))
