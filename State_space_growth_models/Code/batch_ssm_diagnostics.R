## Loops through all models, applies burn-in, and creates final SSM diagnostic figures and parameters.

for(mod_num in 1:19){ ## Change total number of models if needed
  source(here::here('config.R'))
  complete_models <- readRDS(here::here('Outputs', model_list_fn))
  site = complete_models[[mod_num]][['site']]
  pft = complete_models[[mod_num]][['pft']]
  model = complete_models[[mod_num]][['model']]
  chain = complete_models[[mod_num]][['chain']]
  burn_in = complete_models[[mod_num]][['burn_in']]
  
  knitr::spin(here::here('Code/create_ssm_diagnostics.R'))
  
  ## Send HTML files to a common location
  outdir <- here::here('Outputs/Final_SSM_diagnostics')
  file.copy(from = 'create_ssm_diagnostics.html', 
            to = file.path(outdir, 
                           paste0(indir_names[[site]], '_', 
                                  pft, '_', model, '_', chain, '.html')), 
            overwrite = TRUE)
  
  ## Clear memory
  rm(list=ls())
  gc()
}
