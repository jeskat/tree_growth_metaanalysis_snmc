## Data of post-processed statespace model files
date <- '2025-12-09'

## List of meta-regressors
metaregressors <- c('rSDI', #'Heat_load', 'BAperTree_m2ha', 'BA_reduction', 'rSDI_calfire',  
                    'ppt_by_waterYr', "cwd_by_waterYr", "pdsi_by_waterYr", #'tmax_by_waterYr', 
                    "TotalDepth", "AWC_50cm")  # "ksat_50cm","CEC_50cm"

## Mapping of code names to site names
common_nms <- list(
  Blodgett = 'Blodgett',
  LaTour = 'LaTour',
  Sequoia = 'Sequoia',
  STEF = 'STEF',
  Teakettle = 'Teakettle',
  TharpsCreek = "Tharp's Creek",
  WLakeTahoe = 'W. Lake Tahoe',
  UPFUb = 'W. Lake Tahoe',
  Blodgett_FFS15 = 'Blodgett',
  LaTour15 = 'LaTour',
  `STEF-VDT15` = 'STEF',
  Teakettle15 = 'Teakettle',
  Sequoia_FFS15 = 'Sequoia',
  Sequoia_Tharp15 = "Tharp's Creek"
)

## Mapping of PFT names
pft_dict <- list(
  cedar = 'Cedar',
  fir = 'Fir',
  yellow_pine = 'Yellow Pine',
  white_pine = 'Sugar Pine',
  other_conifer = 'Douglas Fir'
)

## Site order (North to South)
site_order <- c("Tharp's Creek", 'Sequoia', 'Teakettle_noOverstory', 
                'Teakettle', 'STEF', 'Blodgett', 
                'W. Lake Tahoe', 'LaTour')


## Fixed effects
fxdEff <- c('int_overall', 'slope_size', 'd2_size', 'betaCWD', 
            'treat_mean[Burn]', 'treat_mean[Thin]', 'treat_mean[Burn+Thin]', 
            'treat_cwd[Burn]', 'treat_cwd[Thin]', 'treat_cwd[Burn+Thin]',
            'treat_size[Burn]', 'treat_size[Thin]', 'treat_size[Burn+Thin]')

## Random effects
othrPrms <- c('int_tree_sd', 'int_plot_sd', 'obs_sd', 'res_sd')

