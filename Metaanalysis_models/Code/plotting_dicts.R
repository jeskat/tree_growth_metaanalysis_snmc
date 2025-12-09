### This file provides mappings for plotting styles and nomenclatures


## PFT colors
pft_cols <- list(
  Cedar = '#117733', 
  `Sugar Pine` = '#CC6677', 
  `Yellow Pine` = '#DDCC77', 
  Fir = '#88CCEE', 
  `Douglas Fir` = '#882255') 

## Treatment linestyles, shapes, and colors
trt_styles <- list(
  None = 1, 
  Burn = 2, 
  Thin = 4, 
  `Burn+Thin` = 3,
  Control = 1)

trt_shapes <- list(
  None = 16, 
  Burn = 17, 
  Thin = 3, 
  `Burn+Thin` = 15, 
  Control = 16)  

trt_cols <- list(None = 'grey', 
                 Burn = '#D55E00', 
                 'Thin' = '#0072B2', 
                 `Burn+Thin` = '#42017C', 
                 Control = '#787276', 
                 `Burn:Thin Interaction` = '#42017C')

## Site shapes
site_shapes <- list(Blodgett = 0, 
                    LaTour = 1, 
                    Sequoia = 2, 
                    STEF = 3, 
                    Teakettle = 4, 
                    `Tharp's Creek` = 5,
                    `W. Lake Tahoe` = 6)

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

## Mapping of parameter names for plots
pretty_nms <- list(
  int_overall = "$\\beta_0$",
  slope_size = "$\\beta_1$",
  d2_size = "$\\beta_2$",
  betaCWD = "$\\beta_3$",
  int_plot_sd = "sd($\\alpha_{unit}$)",
  int_tree_sd = "sd($\\alpha_{tree}$)",
  int_year_sd = "sd($\\alpha_{year}$)",
  int_block_sd = "sd($\\alpha_{block}$)",
  obs_sd = "$\\sigma_{obs}$",
  res_sd = "$\\sigma_{proc}$",
  `treat_mean[BURN]` = "$\\beta_{Burn}$",
  `treat_mean[THIN]` =  "$\\beta_{Thin}$",
  `treat_mean[BURN+THIN]` = "$\\beta_{Burn+Thin}$",
  `treat_mean[Burn]` = "$\\gamma_{Burn}$",
  `treat_mean[Thin]` =  "$\\gamma_{Thin}$",
  `treat_mean[Burn+Thin]` = "$\\gamma_{Burn+Thin}$",
  `treat_cwd[Burn]` = "$\\kappa_{Burn}$",
  `treat_cwd[Thin]` =  "$\\kappa_{Thin}$",
  `treat_cwd[Burn+Thin]` = "$\\kappa_{Burn+Thin}$",
  `treat_size[Burn]` = "$\\nu_{Burn}$",
  `treat_size[Thin]` =  "$\\nu_{Thin}$",
  `treat_size[Burn+Thin]` = "$\\nu_{Burn+Thin}$",
  BTintrctn = "$\\beta_{intrctn}$",
  Avg_thin = "Average effect of thinning",
  Avg_burn = "Average effect of prescribed burning",
  beta3 = "Effect of basal area of trees >= 15cm DBH"
)

## Mapping of metaregressor names
pretty_mrs <- list(
  BA_reduction= 'Treatment intensity (% basal area reduction)',
  Elevation='Elevation (m)',
  tempC= 'Average summer maximum temperature (°C)',
  tmax_by_waterYr = 'Average summer maximum \ntemperature (°C)',
  precipitation_amount= 'Average winter \nprecipitation (mm)',
  ppt_by_waterYr = 'Average winter \nprecipitation (mm)',
  snow_duration_days= 'Average snow duration (days)',
  terraclimate_cwd_annual= 'Average annual CWD (mm)',
  terraclimate_pdsi_annual= 'Average annual PDSI',
  cwd_by_waterYr = 'Average annual CWD (mm)',
  pdsi_by_waterYr = 'Average annual PDSI',
  TotalDepth= 'Soil depth (cm)',
  CEC_50cm= 'CEC in top 50cm (meq/100g)',
  OM_50cm= 'Organic matter \n(% by weight)',
  ksat_50cm= "Hydraulic conductivity \nin top 50 cm (um/s)",
  AWC_50cm = "Available water capacity \nin top 50 cm (cm/cm)",
  AWS150 = "Available water storage in top 150 cm (cm)",
  `Number.of.entries` = 'Number of entries',
 BAperTree_m2ha = 'Mean post-treatment basal area (m2/ha)', 
 Slope = "Slope (%)",
 Folded_aspect = 'Folded aspect (degrees)',
 Folded_northness = 'Folded aspect (degrees)',
 Latitude = 'Latitude (degrees North)',
 Heat_load = 'Heat load index',
 rSDI = 'Relative stand density index',
 rSDI_calfire = 'Relative stand density index'
)

## Site order (North to South)
site_order <- c("Tharp's Creek", 'Sequoia', 'Teakettle_noOverstory', 
                'Teakettle', 'STEF', 'Blodgett', 
                'W. Lake Tahoe', 'LaTour')

## Shorthand metaregressor names for correlation plot
corrplot_shorthand <- list(
  'Site' = 'Site',
  'Treatment' = 'Treatment',
  'BA_reduction' = "Intensity",
  'Elevation' = 'Elevation',
  'tempC' = 'Max summer temp',
  'precipitation_amount' = 'Winter precip',
  'snow_duration_days' = 'Snow duration',
  'terraclimate_cwd_annual' = 'CWD',
  'terraclimate_pdsi_annual' = 'PDSI',
  'tmax_by_waterYr' = 'Max summer temp',
  'ppt_by_waterYr' = 'Winter precip',
  'cwd_by_waterYr' = 'CWD',
  'pdsi_by_waterYr' = 'PDSI',
  'TotalDepth' = 'Soil depth',
  'CEC_50cm' = 'CEC',
  'OM_50cm' = "Organic matter",
  'ksat_50cm' = 'Hydraulic conductivity',
  'AWC_50cm' = 'AWC',
  'AWS150' = 'AWS',
  'Number.of.entries' = 'Number of entries',
  'BAperTree_m2ha' = 'Basal area',
  'Slope' = 'Slope',
  'Folded_aspect' = 'Folded aspect',
  'Folded_northness' = 'Folded aspect',
  'Latitude' = 'Latitude',
  'Heat_load' = 'Heat load',
  'rSDI' = 'rSDI',
  'rSDI_calfire' = 'rSDI'
)

## Order of SSM parameters
ssm_param_order <- c("$\\beta_0$", "$\\beta_1$", "$\\beta_2$", "$\\beta_3$",
                     "$\\gamma_{Burn}$", "$\\gamma_{Thin}$", "$\\gamma_{Burn+Thin}$",
                     "$\\kappa_{Burn}$", "$\\kappa_{Thin}$", "$\\kappa_{Burn+Thin}$",
                     "$\\nu_{Burn}$", "$\\nu_{Thin}$", "$\\nu_{Burn+Thin}$",
                     "sd($\\alpha_{tree}$)", "sd($\\alpha_{unit}$)", 
                     "$\\sigma_{obs}$",  "$\\sigma_{proc}$")

## Fixed effects
fxdEff <- c('int_overall', 'slope_size', 'd2_size', 'betaCWD', 
            'treat_mean[Burn]', 'treat_mean[Thin]', 'treat_mean[Burn+Thin]', 
            'treat_cwd[Burn]', 'treat_cwd[Thin]', 'treat_cwd[Burn+Thin]',
            'treat_size[Burn]', 'treat_size[Thin]', 'treat_size[Burn+Thin]')

## Random effects
othrPrms <- c('int_tree_sd', 'int_plot_sd', 'obs_sd', 'res_sd')

