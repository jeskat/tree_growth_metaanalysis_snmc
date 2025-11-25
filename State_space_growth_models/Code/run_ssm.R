library(stringi)
# install.packages('markdown', repos="http://cran.us.r-project.org")
library(markdown)
library(tictoc)

source(here::here('State_space_growth_models/Code/ssm_config.R'))

###------------------------------------------------------------##
## Get inputs passed from shell script
###------------------------------------------------------------##

args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args[1]))
}

### The following arguments are passed in from args_ssm.sh

### Plant functional type (cedar, fir, yellow_pine, white_pine)
pft <- paste0(args[1]) 
print(paste0("PFT: ", pft))

#### Study site (Blodgett, LaTour, Sequoia, STEF, Teakettle, TharpsCreek, WLakeTahoe)
site <- paste0(args[2]) 
print(paste0('Site: ',site))

### SSM (ddbh_model, ddbh_model_TreePrior, ddbh_model_UnitPrior)
model_to_run <- paste0(args[3]) 
print(paste0('Model: ',model_to_run))

### Suffix for run
note <- args[4]
print(paste0('Note: ',note))

### For two chains: specify c1 or c2
chain <- args[5]
print(paste0('Chain: ', chain))

###------------------------------------------------------------##
## Configure run
###------------------------------------------------------------##

if(chain=='c1'){
  seed = 2024 # Set seed for this chain
    
}else if(chain=='c2'){
  seed = 4202
}
print(paste0('Seed: ', seed))

### Configuration for current run
source(here::here('State_space_growth_models', 'Code', 'scripts_for_parallel_processing', 
                 paste0('config_', chain, '.R'))) 
### List of prior runs in chain
pr_rn_lst <- paste0('runs_to_concatenate_', chain, '.RData') 

### Extract configs for this site and PFT
run_config <- configs[[site]][[pft]]

n_iter <- run_config$n_iter # Number of iterations for MCMC
print(paste0('MCMC iterations: ',n_iter))

thin_interval <- run_config$thin_interval # Thinning interval for MCMC
print(paste0("thin_interval: ",thin_interval))


###------------------------------------------------------------##
## Make an outputs folder for this model run
###------------------------------------------------------------##

# Identify prior runs in this chain
runs <- readRDS(here::here('State_space_growth_models/Outputs', pr_rn_lst))
prior_runs <- runs[[model_to_run]][[site]][[pft]]

# Run immediately prior to the current one, from which to continue
prior_run_nm <- prior_runs[length(prior_runs)]

# Replace file path as needed
output_file_path <- mcmc_results_dir #file.path('Outputs')
prior_run_dir <- file.path(output_file_path, prior_run_nm)
print(paste0('prior_run_dir: ', prior_run_dir))

if(is.null(prior_run_nm) & run_config$first_run==FALSE){
  # Exit if we're looking for a prior run and can't find it.
  print('No existing prior run. Aborting.') 
}else{
  # Create a new directory to store MCMC results from current run.
  new_dir <- stri_replace_all_regex('Z_Y_W_J_K',pattern = c('Z','Y','W','J','K'), 
                                    replacement = c(as.character(Sys.Date()),
                                                    as.character(site), 
                                                    as.character(model_to_run), 
                                                    paste0(n_iter/1000,'k'),
                                                    note),
                                    vectorize=FALSE)
  
  outdir <- file.path(output_file_path,new_dir)
  
  print(paste0('outdir: ', outdir))
  
  if (file.exists(here::here(outdir))){
    print("Overwriting older model of same name")
  } else {
    dir.create(here::here(outdir))
  }
  
  # Output file name
  fileName <- stri_replace_all_regex('Z_Y_X_W_J_K',
                                     pattern = c('Z','Y','X','W','J', 'K'),
                                     replacement = c(as.character(Sys.Date()),
                                                     as.character(site),
                                                     as.character(pft),
                                                     as.character(model_to_run),
                                                     paste0(n_iter/1000,'k'),
                                                     as.character(note)),
                                     vectorize=FALSE)
  
  
  
  ###------------------------------------------------------------##
  ## Process data inputs
  ###------------------------------------------------------------##
  
  #Keep track of the amount of time this run takes. 
  tic()
  
  # Create input dataframes for this site, PFT, and model
  source(here::here('State_space_growth_models', 'Code', 'process_inputs.R'))
  
  
  ###------------------------------------------------------------##
  ## NIMBLE model set-up
  ###------------------------------------------------------------##
  
  # Set random seed for reproducibility 
  set.seed(seed)
  
  # Read nimble model file
  model_file <- stri_replace_all_regex('Z.R', 
                                       pattern = c('Z'), 
                                       replacement = as.character(model_to_run), 
                                       vectorize = FALSE)
  
  source(here::here('State_space_growth_models', 'Code', 'nimble_models', model_file))
  source(here::here('State_space_growth_models', 'Code', 'ssm_restart_functions.R'))
  
  # Create and compile nimble model
  useWAIC <- TRUE
  
  model <- nimbleModel(stateSpaceCode,
                       data = data,
                       constants = constants,
                       inits = inits)
  
  
  mcmcConf <- configureMCMC(model, thin = thin_interval, monitors2 = monitors2, 
                            thin2 = thin_interval, enableWAIC = useWAIC)
  

  mcmc <- buildMCMC(mcmcConf)
  cmodel <- compileNimble(model, showCompilerOutput = TRUE)
  cmcmc <- compileNimble(mcmc, project = model)
  
  ##----------------------------------------------------------------------------
  ## Restart from existing run
  ##----------------------------------------------------------------------------
  if(is.null(prior_run_nm)==FALSE){
    # If MCMC is continuing from a previous run, load state variables
    stateList <- readRDS(file = here::here(prior_run_dir, 
                                          paste0(pft, '_restart_states.RData')))
    
    modelState <- stateList$modelState
    mcmcState <- stateList$mcmcState
    if(useWAIC){
      waicState <- stateList$waicState
    }
    
    ## restore the saved "state" into the new model and new MCMC
    setModelState(cmodel, modelState)
    setMCMCstate(mcmcConf, cmcmc, mcmcState)
    if(useWAIC){
      setWAICstate(cmcmc, waicState)
    }
    
    .Random.seed <- stateList$seed
    
    ## continue MCMC run
    cmcmc$run(niter = n_iter, thin = thin_interval, thin2 = thin_interval,
              reset = FALSE, resetWAIC = FALSE)
    
  }else{
    ## Only for the first run of a chain
    set.seed(seed)
    cmcmc$run(niter = n_iter, thin = thin_interval, thin2 = thin_interval)
  }
  
  
  # ----------------------------------------------------------------------------
  # Save MCMC samples and states such that we can pick up where we left off
  # ----------------------------------------------------------------------------
  
  # For all runs
  samples <- as.matrix(cmcmc$mvSamples)
  samples2 <- as.matrix(cmcmc$mvSamples2)
  
  # Get states at end of current run
  stateList <- list(modelState = getModelState(cmodel),
                    mcmcState = getMCMCstate(mcmcConf, cmcmc),
                    seed = .Random.seed)
  if(useWAIC){
    stateList[['waicState']] <- getWAICstate(cmcmc)}
  
  saveRDS(stateList, file = here::here(outdir, 
                                      paste0(pft, '_restart_states.RData')))
  saveRDS(samples, file = here::here(outdir, paste0(fileName, '.RData')))
  saveRDS(samples2, file = here::here(outdir, paste0(fileName, '_2', '.RData')))
  
  
  print(cmcmc$getWAIC())
  
  # End timer
  toc()
  
  # Move log file to the appropriate directory
  log_fn <- paste0(pft, '_', site, '_', model_to_run, '_', note, '.Rout')
  
  if(file.copy(from = log_fn, to = here::here(outdir, log_fn), overwrite = TRUE)) {
    # Remove file if copy is successful
    file.remove(log_fn)
    print("File moved successfully using copy/delete method.")
  }else{
    print('Log file copy failed')
  }
  
  # Concatenate new run to the list of previous runs
  runs_current <- readRDS(here::here('State_space_growth_models/Outputs', pr_rn_lst))
  
  runs_current[[model_to_run]][[site]][[pft]] <- append(runs[[model_to_run]][[site]][[pft]], new_dir)
  saveRDS(runs_current, here::here('State_space_growth_models/Outputs', pr_rn_lst))
  
  # Create results summary and diagnostics
  model <- model_to_run
  knitr::spin(here::here('State_space_growth_models', 'Code', 'create_ssm_diagnostics.R'))
  
  # Save in appropriate directory
  if(file.copy(from = 'create_ssm_diagnostics.html', 
               to = here::here(outdir, 
                               paste0(site, '_', pft, '_', model, '_', chain, '.html')), 
               overwrite = TRUE)) {
    # Remove file if copy is successful
    file.remove('create_ssm_diagnostics.html')
    print("File moved successfully using copy/delete method.")
  }else{
    print('Diagnostics file copy failed')
  }
}