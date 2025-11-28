###----------------------------------------------------------------------------##
## Loop through all models and compile Gelman-Rubin statistics for all parameters
###----------------------------------------------------------------------------##

source(here::here('config.R'))
complete_models <- readRDS(here::here('Outputs', model_list_fn))

## Instantiate empty list
grdf <- vector('list', length = length(complete_models))

for(mod_num in 1:length(complete_models)){ 
  site <- complete_models[[mod_num]][['site']] ## This needs to be mod_site, not site
  model <- complete_models[[mod_num]][['model']]
  pft <- complete_models[[mod_num]][['pft']]
  burn_in <- complete_models[[mod_num]][['burn_in']]
  
  ## Get Gelman-Rubin diagnostics
  source(here::here('Code/convergence_diagnostics_one_model.R'))
  
  ## Dataframe of point estimates
  to_bind <- g$psrf
  to_bind <- data.frame(to_bind)
  to_bind$site <- site
  to_bind$pft <- pft
  to_bind$model <- model
  to_bind$param <- rownames(to_bind)
  rownames(to_bind) <- NULL
  
  grdf[[mod_num]] <- to_bind
  
}


## List to df 
out <- do.call(rbind, grdf)

## Save
write.csv(out, here::here('Outputs', 
                          paste0('all_gr_diagnostics_', Sys.Date(), '.csv')), 
          row.names = FALSE)
