### Randomly sample three trees from a data frame of tagged individuals
sample_trees <- function(df, size){
  mini_df <- df[df$Size_class==size,]
  samps <- sample(1:dim(mini_df)[1],3)
  return(rownames(mini_df[samps,]))
}

### Function to identify which column name includes the log(DBH) parameter.
## Need this because some of the SSM models in Katz et al. were run using
## "slope_size" as the name of the log(DBH) parameter, rather than "log_size."
get_lnD_colnm <- function(mcmc_outputs){
  if('slope_size' %in% colnames(mcmc_outputs)){
    lnD_nm <- 'slope_size'
  }else if('log_size' %in% colnames(mcmc_outputs)){
    lnD_nm <- 'log_size'
  }else{
    print('log(DBH) parameter is missing from results.')
  }
  return(lnD_nm)
}


### Load and combine into one data table all MCMC iterations associated with a chain
### mcmc_results_dir is where the subdirectories holding the MCMC outputs are saved
load_mcmc <- function(site, pft, model, burn_in, chain, mcmc_results_dir, to_report){
  ## Read in the list of prior runs in this chain.
  if(chain == 'c1'){
    runs <- readRDS(here::here('State_space_growth_models/Outputs', runs_c1_fn))
  }else if(chain == 'c2'){
    runs <- readRDS(here::here('State_space_growth_models/Outputs', runs_c2_fn))
  }
  
  ## List, in chronological order, of the output sub-directories that contain the 
  ## MCMC iterations for this site and PFT.
  outdirs <- runs[[model]][[site]][[pft]]
  
  ## Load and combine into one data table all MCMC iterations associated with this chain
  for(i in 1:length(outdirs)){ #Iterate through runs in chain
    # Get the name of the .RData file associated with the run
    substrs <- unlist(strsplit(outdirs[i], '_'))
    fn <- paste(c(substrs[1], site, pft, model, substrs[(length(substrs)-1)], 
                  substrs[length(substrs)]), collapse = '_')
    
    if(file.exists(here::here(mcmc_results_dir, outdirs[i], paste0(fn, '.RData')))==FALSE){
      fn <- paste(c(substrs[1], site, pft, model, substrs[(length(substrs)-1)]), 
                  collapse = '_')}
    
    # Open .RData file and rbind to previous run
    if(i==1){ # For first run only
      all_samples <- readRDS(here::here(mcmc_results_dir, outdirs[i], paste0(fn, '.RData')))
      print(paste0('Dimensions of run ', i, ": ", dim(all_samples)))
      print(all_samples[(nrow(all_samples)-3):nrow(all_samples), 1:10])
      
      # Identify which column name includes the log(DBH) parameter
      lnD_nm <- get_lnD_colnm(all_samples)
      
      # add parameters to report
      to_report <- c(to_report, 
                     lnD_nm,
                     colnames(all_samples)[grep('beta', colnames(all_samples))], # all betas
                     colnames(all_samples)[grep('_sd', colnames(all_samples))], # all SDs
                     colnames(all_samples)[grep('treat_',colnames(all_samples))] # all treatment effects and interactions
      )
      
      #Subset columns to report
      if(all(is.na(all_samples[1,]))){
        print('First record is null, what happened???')
        all_samples <- all_samples[2:dim(all_samples)[1], to_report]
      }else{
        all_samples <- all_samples[, to_report]
      }
      
      
    }else{ # for all subsequent runs
      new_samples <- readRDS(here::here(mcmc_results_dir, outdirs[i], paste0(fn, '.RData')))
      print(paste0('Dimensions of run ', i, ": ", dim(new_samples)))
      ## Print first few row of samples (select vars)
      print(new_samples[1:3, 1:10])
      
      # Discard the first row if it contains only NA
      if(all(is.na(new_samples[1,]))){
        all_samples <- rbind(all_samples, new_samples[2:dim(new_samples)[1], to_report])
      }else{
        print('Warning: connected run did not start with NA. Should this be the first run?')
        all_samples <- rbind(all_samples, new_samples[, to_report])
      } 
    }
    
  }
  
  # Rename treat_means 
  treat_cols <- to_report[grep('treat_mean', to_report)]
  
  for(i in parse_number(treat_cols)){
    for(n in c('treat_mean[', 'treat_cwd[', 'treat_size[')){
      old_nm <- paste0(n, i, ']')
      trt_nm <- levels(plotData$fTreatment)[i]
      new_nm <- paste0(n, trt_nm, ']')
      colnames(all_samples)[colnames(all_samples) == old_nm] <- new_nm
      to_report[to_report == old_nm] <- new_nm
    }
  }
  
  ## How many samples total?
  nSamps <- nrow(all_samples)
  
  ## Apply burn-in
  ## Create an MCMC object for coda functions
  all_samples <- as.mcmc(as.matrix(all_samples[burn_in:nSamps,]))
  
  return(all_samples)
}

### Load and combine into one data table all MCMC iterations associated with a chain
### Similar to load_mcmc, but specifically gets plot random effects
load_random_plot_effects <- function(site, pft, model, burn_in, chain, mcmc_results_dir){
  
  ## Read in the list of prior runs in this chain.
  if(chain == 'c1'){
    runs <- readRDS(here::here('State_space_growth_models/Outputs', runs_c1_fn))
  }else if(chain == 'c2'){
    runs <- readRDS(here::here('State_space_growth_models/Outputs', runs_c2_fn))
  }
  
  ## List, in chronological order, of the output sub-directories that contain the 
  ## MCMC iterations for this site and PFT.
  outdirs <- runs[[model]][[site]][[pft]]  
  
  ## Load and combine into one data table all MCMC iterations associated with this chain
  for(i in 1:length(outdirs)){ #Iterate through runs in chain
    # Get the name of the .RData file associated with the run
    substrs <- unlist(strsplit(outdirs[i], '_'))
    fn <- paste(c(substrs[1], site, pft, model, substrs[(length(substrs)-1)], 
                  substrs[length(substrs)]), collapse = '_')
    
    if(file.exists(here::here(mcmc_results_dir, outdirs[i], paste0(fn, '.RData')))==FALSE){
      fn <- paste(c(substrs[1], site, pft, model, substrs[(length(substrs)-1)]), 
                  collapse = '_')}
    
    ## Get random effects
    # Open .RData file and rbind to previous run
    if(i==1){ # For first run only
      random_effect_samples <- readRDS(here::here(mcmc_results_dir, outdirs[i], paste0(fn, '_2.RData')))
            
      #Subset columns to report
      nms <- colnames(random_effect_samples)[grep('int_plot', colnames(random_effect_samples))]
      random_effect_samples <- random_effect_samples[, nms]
      
      
    }else{ # for all subsequent runs
      new_rndm_fx <- readRDS(here::here(mcmc_results_dir, outdirs[i], paste0(fn, '_2.RData')))
      
      # Discard the first row if it contains only NA
      if(all(is.na(new_rndm_fx[1,]))){
        random_effect_samples <- rbind(random_effect_samples, new_rndm_fx[2:dim(new_rndm_fx)[1],nms])
      }else{
        print('Warning: connected run did not start with NA. Should this be the first run?')
        random_effect_samples <- rbind(random_effect_samples, new_rndm_fx[,nms])
        } 
      }
    }
  
  ## Apply burn-in, remove columns with all zeros
  plot_effects <- random_effect_samples[burn_in:nrow(random_effect_samples), nms]
  clean_plot_effects <- plot_effects[, which(colSums(plot_effects)!=0)]
  
  
  return(clean_plot_effects)
}


### Report summary statistics for the posterior MCMC of a given parameter
summarize_mcmc <- function(out, param){
  tryCatch(
    {
      out_met <- summary(out[,param])
      avg <- out_met$statistics['Mean']
      low_ci <- out_met$quantiles['2.5%']
      up_ci <- out_met$quantiles['97.5%']
      stdv <- out_met$statistics['SD']
      naiveSE <- out_met$statistics['Naive SE']
      tsSE <- out_met$statistics['Time-series SE']
      return(c(param, avg, low_ci, up_ci, stdv, naiveSE, tsSE))
    }, 
    error=function(e){
      message(paste('An error occurred:', param))
      print(e)
    },
    warning=function(w){
      message(paste('A warning occurred:', param))
      print(w)
      return(NA)
    }
  )
}


## Function to calculate the mean values of a set of posterior distributions
### Takes in parameter names (nms) and an MCMC table-like object in which 
### each column contains a with posteriors associated with one of the parameters.
summarize_effect_by_trt <- function(nms, mcmc_by_trt){
  ests <- numeric(length(nms))
  
  for(i in 1:length(nms)){
    var <- summarize_mcmc(mcmc_by_trt, nms[i])
    ests[i] <- as.numeric(var['Mean'])
  }
  
  return(data.frame(variable=nms, mean=ests))
}



### Function to add columns for site, PFT, and model
addSitePFTModCols <- function(d, site, pft, model){
  d$site <- rep(site, nrow(d))
  d$pft <- rep(pft, nrow(d))
  d$model <- rep(model, nrow(d))
  return(d)
}

get_reparameterized_plot_growth <- function(CWD, # Climatic water deficit in units of standard deviations from site mean
                                           DBH, # DBH in units of cm
                                           plot_ids, # 
                                           plotData, 
                                           all_samples,
                                           obs_ln_mean,
                                           obs_ln_stdv,
                                           obs_sq_mean,
                                           obs_sq_sd){
  ## Instantiate lists to hold MCMC summary statistics for each plot
  out <- vector('list', length = 0)
  plts <- vector('numeric', length = length(plot_ids))
  means <- vector('numeric', length = length(plot_ids))
  sds <- vector('numeric', length = length(plot_ids))
  treats <- vector('character', length(plot_ids))
  
  # Identify which column name includes the log(DBH) parameter
  lnD_nm <- get_lnD_colnm(all_samples)
  
  for(p in 1:length(plot_ids)){
    # Get treatment associated with plot
    trt <- plotData[as.character(plotData$PlotID) == plot_ids[p], 'treatment']
    nm <- paste0('treat_mean[', trt, ']')
    treats[p] <- nm
    plts[p] <- p
    
    ## Calculate growth in plot at the specified DBH and CWD
    out[[p]] <- all_samples[, 'int_overall'] + ## overall intercept
      all_samples[, nm] + ## treatment effect
      all_samples[, lnD_nm] * (log(DBH) - obs_ln_mean)/obs_ln_stdv + ## tree size effects
      all_samples[, 'd2_size'] * (DBH^2 - obs_sq_mean)/obs_sq_sd + 
      all_samples[, 'betaCWD']*CWD + ## CWD effect
      all_samples[, paste0('treat_cwd[', trt, ']')] * CWD + ## CWD:treatment interaction
      all_samples[, paste0('treat_size[', trt, ']')] * (log(DBH) - obs_ln_mean)/obs_ln_stdv +
      plot_effects[, paste0('int_plot[', p, ']')] ## random plot effect
    
    ## Save MCMC means and SDs to lists
    means[p] <- summary(out[[p]])$statistics['Mean']
    sds[p] <- summary(out[[p]])$statistics['SD']
  }
  
  df <- data.frame(out)
  names(df) <- plot_ids
  
  long_tbl <- data.frame(unit =plot_ids, 
                         variable = treats, 
                         mean=means, 
                         sd = sds,
                         DBH = DBH,
                         CWD = CWD)
  
  ## Summary table with plot ID, treatment, and mean and SD of growth
  return(list(long_tbl = long_tbl, cv_mat = cov(df))) 
}

