###---------------------------------------------------------------------------##
## Runs meta-regressions for climate, soil, and stand structure variables  
###---------------------------------------------------------------------------##


# Dependencies ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggpattern)
library(latex2exp)
library(tidyverse)
library(metafor)
library(Matrix)
library(clubSandwich)
library(scico)

source(here::here('Metaanalysis_models/Code/metaanalysis_config.R'))
source(here::here('Metaanalysis_models/Code/metaanalysis_functions.R'))

source(here::here('Figures_and_Tables/Code/plotting_config.R'))
source(here::here('Figures_and_Tables/Code/plotting_functions.R'))

### Load meta-regressors
rgrs_by_unit_raw <- read.csv(
  here::here('Covariate_and_metaregressor_processing/Processed_data/all_regressors_by_unit.csv'))

## Load growth effects by unit (reparameterized to DBH = 30cm, CWD = 0)
growth_effects <- read.csv(here::here('Metaanalysis_models/Reparameterized_SSM_outputs/', 
                                      paste0('unit_growth_effects_by_trt_',
                                             date, 
                                             '.csv')))

## Load covariance matrices associated with fixed_effects (unit-level)
growth_cov_unit <- readRDS(here::here('Metaanalysis_models/Reparameterized_SSM_outputs/',
                                      paste0('cov_unit_growth_by_trt_', 
                                             date, 
                                             '.RData')))

# Format data -------------------------------------------------------------

growth_effects$Treatment <- sub(".*\\[([^][]+)].*", "\\1", 
                                growth_effects$variable)
unit_effects <- formatForMetafor(growth_effects)
unit_cov_m <- lapply(growth_cov_unit, function(x) as.matrix(x))

## Format site names
rgrs_by_unit_raw$common_nm <- NaN

for(i in 1:length(rgrs_by_unit_raw$Site)){
  rgrs_by_unit_raw$common_nm[i] <- common_nms[rgrs_by_unit_raw$Site[i]]
}

rgrs_by_unit_raw$common_nm <- factor(rgrs_by_unit_raw$common_nm, site_order)



# Relative stand density index (rSDI) -------------------------------------

mr_rsdi <- make_meta_df(full_data = unit_effects, 
                        cv_list = unit_cov_m,
                        formula = ~(Burn*Thin*rSDI), 
                        by_unit = TRUE,
                        mtrgrsn = TRUE,
                        regressor = 'rSDI',
                        rgr_df = rgrs_by_unit_raw)


rsdi_tbl <- make_results_table(mr_rsdi[[2]], by_unit = TRUE)
rsdi_tbl <- beautify_table_labels(rsdi_tbl, 'rSDI')

rsdi_plot <- plotMetaRegression(mr_rsdi,
                                rgr = 'rSDI', 
                                sgnf = sgnf, 
                                by_unit=TRUE, 
                                uncenter = TRUE) 
print(rsdi_plot)

write.csv(rsdi_tbl, here::here(tbl_dir, 'mtrgr_rsdi.csv'), row.names = FALSE)

ggsave(here::here(fig_dir, 'mtrgr_rsdi.jpeg'),
       rsdi_plot +coord_cartesian(clip = 'off'), device = 'jpeg',
       dpi = 'print',
       width = 6.5, height = 8, units = 'in')


# rSDI without treatment effects------------------------------------------------

## rSDI regression without a treatment variable
mr_rsdi_noTrt <- make_meta_df(full_data = unit_effects, 
                                         cv_list = unit_cov_m,
                                         formula = ~(rSDI), 
                                         by_unit = TRUE,
                                         mtrgrsn = TRUE,
                                         regressor = 'rSDI',
                                         rgr_df = rgrs_by_unit_raw)


mr_rsdi_noTrt_tbl <- make_results_table(mr_rsdi_noTrt[[2]], by_unit = TRUE)
mr_rsdi_noTrt_tbl <- beautify_table_labels(mr_rsdi_noTrt_tbl, 'rSDI')

write.csv(mr_rsdi_noTrt_tbl, here::here(tbl_dir, 'mtrgr_rsdi_noTrt.csv'), 
          row.names = FALSE)

# rSDI without a treatment:rSDI interaction-------------------------------------

mr_rsdi_noIrctn <- make_meta_df(full_data = unit_effects, 
                        cv_list = unit_cov_m,
                        formula = ~(Burn*Thin+rSDI), 
                        by_unit = TRUE,
                        mtrgrsn = TRUE,
                        regressor = 'rSDI',
                        rgr_df = rgrs_by_unit_raw)


rsdi_noIrctn_tbl <- make_results_table(mr_rsdi_noIrctn[[2]], by_unit = TRUE)
rsdi_noIrctn_tbl <- beautify_table_labels(rsdi_noIrctn_tbl, 'rSDI')

write.csv(rsdi_noIrctn_tbl, here::here(tbl_dir, 'mtrgr_rsdi_noIrctn.csv'), 
          row.names = FALSE)



# rSDI (Burn only) --------------------------------------------------------

## Remove all thinning treatments and sites that did not burn
burn_only <- unit_effects[unit_effects$Treatment %in% c('None', 'Burn') & 
                      !(unit_effects$common_nm %in% c('W. Lake Tahoe', 'LaTour')),]

burn_units <- unique(burn_only$unit)

cov_list_burn <- lapply(unit_cov_m[!(grepl('UPFUb', names(unit_cov_m))) & 
                                     !(grepl('LaTour', names(unit_cov_m)))],
                        function(x){nms <- which(colnames(x) %in% burn_units)
                        x[nms,nms]})

mr_rsdi_burn <- make_meta_df(full_data = burn_only, 
                                cv_list = cov_list_burn,
                                formula = ~(Burn*rSDI), 
                                by_unit = TRUE,
                                mtrgrsn = TRUE,
                                regressor = 'rSDI',
                                rgr_df = rgrs_by_unit_raw)


rsdi_burnOnly_tbl <- make_results_table(mr_rsdi_burn[[2]], by_unit = TRUE)
rsdi_burnOnly_tbl <- beautify_table_labels(rsdi_burnOnly_tbl, 'rSDI')

write.csv(rsdi_burnOnly_tbl, here::here(tbl_dir, 'mtrgr_rsdi_burnOnly.csv'), 
          row.names = FALSE)



# rSDI (Thin only) --------------------------------------------------------

## Remove all burn treatments and sites that did not thin
thin_only <- unit_effects[unit_effects$Treatment %in% c('None', 'Thin') & 
                      !(unit_effects$common_nm %in% c('Sequoia', "Tharp's Creek")),]

thin_units <- unique(thin_only$unit)

cov_list_thin <- lapply(unit_cov_m[!(grepl('Sequoia_FFS', names(unit_cov_m))) &
                                     !(grepl('Sequoia_Tharp', names(unit_cov_m)))],
                        function(x){nms <- which(colnames(x) %in% thin_units)
                        x[nms,nms]})

## Adjust the data above for STEF (due to overlap in unit names for STEF and LaTour)
thin_site <- unit_effects[unit_effects$Treatment %in% c('None', 'Thin') & 
                            unit_effects$common_nm == 'STEF',]
thin_units_STEF <- unique(thin_site$unit)
rplcmnt <- unit_cov_m[grepl('STEF', names(unit_cov_m))]
rplcmnt <- lapply(rplcmnt, function(x){nms <- which(colnames(x) %in% thin_units_STEF)
x[nms,nms]})

for(n in names(rplcmnt)){
  cov_list_thin[[n]] = rplcmnt[[n]]
}


mr_rsdi_thin <- make_meta_df(full_data = thin_only, 
                             cv_list = cov_list_thin,
                             formula = ~(Thin*rSDI), 
                             by_unit = TRUE,
                             mtrgrsn = TRUE,
                             regressor = 'rSDI',
                             rgr_df = rgrs_by_unit_raw)


rsdi_thinOnly_tbl <- make_results_table(mr_rsdi_thin[[2]], by_unit = TRUE)
rsdi_thinOnly_tbl <- beautify_table_labels(rsdi_thinOnly_tbl, 'rSDI')

write.csv(rsdi_thinOnly_tbl, here::here(tbl_dir, 'mtrgr_rsdi_thinOnly.csv'), 
          row.names = FALSE)



# Winter precipitation ----------------------------------------------------

mr_ppt <- make_meta_df(full_data = unit_effects, 
                       cv_list = unit_cov_m,
                       formula = ~(Burn*Thin*ppt_by_waterYr), 
                       by_unit = TRUE,
                       mtrgrsn = TRUE,
                       regressor = 'ppt_by_waterYr',
                       rgr_df = rgrs_by_unit_raw)


ppt_tbl <- make_results_table(mr_ppt[[2]], by_unit = TRUE)
ppt_tbl <- beautify_table_labels(ppt_tbl, 'ppt_by_waterYr')

ppt_plot <- plotMetaRegression(mr_ppt,
                               rgr = 'ppt_by_waterYr', 
                               sgnf = sgnf, 
                               by_unit=TRUE, 
                               uncenter = TRUE) 
print(ppt_plot)

write.csv(ppt_tbl, here::here(tbl_dir, 'mtrgr_ppt.csv'), row.names = FALSE)

ggsave(here::here(fig_dir, 'mtrgr_ppt.jpeg'),
       ppt_plot +coord_cartesian(clip = 'off') + 
         xlab('Average winter precipitation (mm)'), 
       device = 'jpeg', dpi = 'print',
       width = width, height = height, scale = scale, units = 'in')


# Climate Water Deficit---------------------------------------------------------

mr_cwd <- make_meta_df(full_data = unit_effects, 
                       cv_list = unit_cov_m,
                       formula = ~(Burn*Thin*cwd_by_waterYr), 
                       by_unit = TRUE,
                       mtrgrsn = TRUE,
                       regressor = 'cwd_by_waterYr',
                       rgr_df = rgrs_by_unit_raw)


cwd_tbl <- make_results_table(mr_cwd[[2]], by_unit = TRUE)
cwd_tbl <- beautify_table_labels(cwd_tbl, 'cwd_by_waterYr')

cwd_plot <- plotMetaRegression(mr_cwd,
                               rgr = 'cwd_by_waterYr', 
                               sgnf = sgnf, 
                               by_unit=TRUE, 
                               uncenter = TRUE) 
print(cwd_plot)

write.csv(cwd_tbl, here::here(tbl_dir, 'mtrgr_cwd.csv'), row.names = FALSE)

ggsave(here::here(fig_dir, 'mtrgr_cwd.jpeg'),
       cwd_plot +coord_cartesian(clip = 'off'), 
       device = 'jpeg', dpi = 'print',
       width = width, height = height, scale = scale, units = 'in')


# Palmer Drought Severity Index-------------------------------------------------

mr_pdsi <- make_meta_df(full_data = unit_effects, 
                        cv_list = unit_cov_m,
                        formula = ~(Burn*Thin*pdsi_by_waterYr), 
                        by_unit = TRUE,
                        mtrgrsn = TRUE,
                        regressor = 'pdsi_by_waterYr',
                        rgr_df = rgrs_by_unit_raw)


pdsi_tbl <- make_results_table(mr_pdsi[[2]], by_unit = TRUE)
pdsi_tbl <- beautify_table_labels(pdsi_tbl, 'pdsi_by_waterYr')

pdsi_plot <- plotMetaRegression(mr_pdsi,
                                rgr = 'pdsi_by_waterYr', 
                                sgnf = sgnf, 
                                by_unit=TRUE, 
                                uncenter = TRUE) 
print(pdsi_plot)

write.csv(pdsi_tbl, here::here(tbl_dir, 'mtrgr_pdsi.csv'), row.names = FALSE)

ggsave(here::here(fig_dir, 'mtrgr_pdsi.jpeg'),
       pdsi_plot +coord_cartesian(clip = 'off'), 
       device = 'jpeg', dpi = 'print',
       width = width, height = height, scale = scale, units = 'in')


# Soil depth--------------------------------------------------------------------

mr_depth <- make_meta_df(full_data = unit_effects, 
                         cv_list = unit_cov_m,
                         formula = ~(Burn*Thin*TotalDepth), 
                         by_unit = TRUE,
                         mtrgrsn = TRUE,
                         regressor = 'TotalDepth',
                         rgr_df = rgrs_by_unit_raw)


depth_tbl <- make_results_table(mr_depth[[2]], by_unit = TRUE)
depth_tbl <- beautify_table_labels(depth_tbl, 'TotalDepth')

depth_plot <- plotMetaRegression(mr_depth,
                                 rgr = 'TotalDepth', 
                                 sgnf = sgnf, 
                                 by_unit=TRUE, 
                                 uncenter = TRUE) 
print(depth_plot)

write.csv(depth_tbl, here::here(tbl_dir, 'mtrgr_depth.csv'), row.names = FALSE)

ggsave(here::here(fig_dir, 'mtrgr_depth.jpeg'),
       depth_plot +coord_cartesian(clip = 'off'), 
       device = 'jpeg', dpi = 'print',
       width = width, height = height, scale = scale, units = 'in')

# Available water capacity (AWC) in the top 50cm of soil------------------------

mr_awc <- make_meta_df(full_data = unit_effects, 
                       cv_list = unit_cov_m,
                       formula = ~(Burn*Thin*AWC_50cm), 
                       by_unit = TRUE,
                       mtrgrsn = TRUE,
                       regressor = 'AWC_50cm',
                       rgr_df = rgrs_by_unit_raw)


awc_tbl <- make_results_table(mr_awc[[2]], by_unit = TRUE)
awc_tbl <- beautify_table_labels(awc_tbl, 'AWC_50cm')

awc_plot <- plotMetaRegression(mr_awc,
                               rgr = 'AWC_50cm', 
                               sgnf = sgnf, 
                               by_unit=TRUE, 
                               uncenter = TRUE) 
print(awc_plot)

write.csv(awc_tbl, here::here(tbl_dir, 'mtrgr_awc.csv'), row.names = FALSE)

ggsave(here::here(fig_dir, 'mtrgr_awc.jpeg'),
       awc_plot +coord_cartesian(clip = 'off') + 
         xlab('Available water capacity in top 50 cm (cm/cm)'), 
       device = 'jpeg', dpi = 'print',
       width = width, height = height, scale = scale, units = 'in')
