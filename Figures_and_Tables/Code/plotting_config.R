### This file provides mappings for plotting styles and nomenclatures

## File directories for figures and tables
fig_dir <- 'Figures_and_Tables/Figures/'
tbl_dir <- 'Figures_and_Tables/Tables/'

## Width and height of metaregression plots
width = 6.5
height = 4
scale = 1.25

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
  `Thin+Burn` = 3,
  Control = 1)

trt_shapes <- list(
  None = 16, 
  Burn = 17, 
  Thin = 3, 
  `Burn+Thin` = 15, 
  `Thin+Burn` = 15,
  Control = 16)  

trt_cols <- list(None = 'grey', 
                 Burn = '#D55E00', 
                 'Thin' = '#0072B2', 
                 `Burn+Thin` = '#42017C',
                 `Thin+Burn` = '#42017C',
                 Control = '#787276', 
                 `Burn:Thin Interaction` = '#42017C',
                 `Thin:Burn Interaction` = '#42017C',
                 `Thin x Burn` = '#42017C')


## Site shapes
site_shapes <- list(Blodgett = 0, 
                    LaTour = 1, 
                    Sequoia = 2, 
                    STEF = 3, 
                    Teakettle = 4, 
                    `Tharp's Creek` = 5,
                    `W. Lake Tahoe` = 6)

## Mapping of PFT names
pft_dict <- list(
  cedar = 'Cedar',
  fir = 'Fir',
  yellow_pine = 'Yellow Pine',
  white_pine = 'Sugar Pine',
  other_conifer = 'Douglas Fir'
)

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

## Site order (North to South)
site_order <- c("Tharp's Creek", 'Sequoia', 
                'Teakettle', 'STEF', 'Blodgett', 
                'W. Lake Tahoe', 'LaTour')


## Fixed effects
fxdEff <- c('int_overall', 'slope_size', 'log_size', 'd2_size', 'betaCWD', 
            'treat_mean[Burn]', 'treat_mean[Thin]', 'treat_mean[Burn+Thin]', 
            'treat_cwd[Burn]', 'treat_cwd[Thin]', 'treat_cwd[Burn+Thin]',
            'treat_size[Burn]', 'treat_size[Thin]', 'treat_size[Burn+Thin]')

## Random effects
othrPrms <- c('int_tree_sd', 'int_plot_sd', 'obs_sd', 'res_sd')


## Mapping of parameter names for plots
pretty_nms <- list(
  int_overall = "$\\beta_0$",
  slope_size = "$\\beta_1$",
  log_size = "$\\beta_1$",
  ln_size = "$\\beta_1$",
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
  `treat_mean[BURN+THIN]` = "$\\beta_{Thin+Burn}$",
  `treat_mean[Burn]` = "$\\gamma_{Burn}$",
  `treat_mean[Thin]` =  "$\\gamma_{Thin}$",
  `treat_mean[Burn+Thin]` = "$\\gamma_{Thin+Burn}$",
  `treat_cwd[Burn]` = "$\\kappa_{Burn}$",
  `treat_cwd[Thin]` =  "$\\kappa_{Thin}$",
  `treat_cwd[Burn+Thin]` = "$\\kappa_{Thin+Burn}$",
  `treat_size[Burn]` = "$\\nu_{Burn}$",
  `treat_size[Thin]` =  "$\\nu_{Thin}$",
  `treat_size[Burn+Thin]` = "$\\nu_{Thin+Burn}$",
  BTintrctn = "$\\beta_{intrctn}$",
  Avg_thin = "Average effect of thinning",
  Avg_burn = "Average effect of prescribed burning",
  beta3 = "Effect of basal area of trees >= 15cm DBH"
)

## Mapping of metaregressor names
pretty_mrs <- list(
  BA_reduction= 'Treatment intensity (% basal area reduction)',
  tmax_by_waterYr = 'Average summer maximum \ntemperature (°C)',
  ppt_by_waterYr = 'Average winter \nprecipitation (mm)',
  cwd_by_waterYr = 'Average annual CWD (mm)',
  pdsi_by_waterYr = 'Average annual PDSI',
  TotalDepth= 'Soil depth (cm)',
  AWC_50cm = "Available water capacity \nin top 50 cm (cm/cm)",
  rSDI = 'Relative stand density index'
)

## Shorthand metaregressor names for correlation plot
corrplot_shorthand <- list(
  'Site' = 'Site',
  'Treatment' = 'Treatment',
  'BA_reduction' = "Intensity",
  'tmax_by_waterYr' = 'Max summer temp',
  'ppt_by_waterYr' = 'Winter precip',
  'cwd_by_waterYr' = 'CWD',
  'pdsi_by_waterYr' = 'PDSI',
  'TotalDepth' = 'Soil depth',
  'AWC_50cm' = 'AWC',
  'rSDI' = 'rSDI'
)

## Order of SSM parameters
ssm_param_order <- c("$\\beta_0$", "$\\beta_1$", "$\\beta_2$", "$\\beta_3$",
                     "$\\gamma_{Burn}$", "$\\gamma_{Thin}$", "$\\gamma_{Thin+Burn}$",
                     "$\\kappa_{Burn}$", "$\\kappa_{Thin}$", "$\\kappa_{Thin+Burn}$",
                     "$\\nu_{Burn}$", "$\\nu_{Thin}$", "$\\nu_{Thin+Burn}$",
                     "sd($\\alpha_{tree}$)", "sd($\\alpha_{unit}$)", 
                     "$\\sigma_{obs}$",  "$\\sigma_{proc}$")
