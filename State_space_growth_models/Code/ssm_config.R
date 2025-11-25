### File paths

## Directory containing all directories of MCMC runs
mcmc_results_dir = file.path('/global','scratch','users','jessicakatz','forest_growth_synthesis','Outputs', 'SS_outputs') #file.path('/global','scratch','users','jessicakatz','forest_growth_synthesis','Outputs','Final_SSM_outputs')

## List of sub-directories associated with each model
runs_c1_fn <- 'runs_to_concatenate_c1.RData' #'authors_runs_to_concatenate_c1.RData' # Chain 1
runs_c2_fn <- 'runs_to_concatenate_c2.RData' #'authors_runs_to_concatenate_c2.RData' # Chain 2

## List of complete models
model_list_fn <- 'authors_complete_models.RData'

## log(DBH) effect parameter name in SSM
lnD_nm <- 'log_size' # slope_size if processing author's results; 'log_size' if running a clean set of SSMs

# Base list the parameters to report out on
to_report <- c(
  'int_overall', 
  'd2_size',
  lnD_nm 
)

## Maps different aliases for sites to the appropriate input data folder
indir_names <- list(
  Blodgett_FFS = 'Blodgett',
  LaTour = 'LaTour',
  UPFUb = 'WLakeTahoe',
  `STEF-VDT` = 'STEF',
  Teakettle = 'Teakettle',
  Sequoia_FFS = 'Sequoia',
  Sequoia_Tharp = "TharpsCreek",
  Blodgett_FFS15 = 'Blodgett',
  LaTour15 = 'LaTour',
  `STEF-VDT15` = 'STEF',
  Teakettle15 = 'Teakettle',
  Sequoia_FFS15 = 'Sequoia',
  Sequoia_Tharp15 = "TharpsCreek",
  Blodgett = 'Blodgett',
  LaTour = 'LaTour',
  WLakeTahoe = 'WLakeTahoe',
  STEF = 'STEF',
  Teakettle = 'Teakettle',
  Sequoia = 'Sequoia',
  TharpsCreek = "TharpsCreek"
)
