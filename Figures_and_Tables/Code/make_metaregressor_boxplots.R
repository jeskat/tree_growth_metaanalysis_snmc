###---------------------------------------------------------------------------##
## Makes supplementary figures, including boxplots of unit-level soil and
## stand structure covariates, and scatterplots of site-level climate variables
###---------------------------------------------------------------------------##


# Data and dependencies ---------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

source(here::here('Figures_and_Tables/Code/plotting_config.R'))
source(here::here('Figures_and_Tables/Code/plotting_functions.R'))

### Load meta-regressors
rgrs_by_unit_raw <- read.csv(
  here::here('Covariate_and_metaregressor_processing/Processed_data/all_regressors_by_unit.csv'))


# Format ------------------------------------------------------------------

## Rename meta-regressor columns
new_colnames <- vector('character', length(colnames(rgrs_by_unit_raw)))
for(i in 1:length(colnames(rgrs_by_unit_raw))){
  if(colnames(rgrs_by_unit_raw)[i] %in% names(pretty_mrs)){
    new_colnames[i] <- pretty_mrs[[colnames(rgrs_by_unit_raw)[i]]]
  }else{
    new_colnames[i] <- colnames(rgrs_by_unit_raw)[i]
  }
}

rgrs_for_splmnt <- rgrs_by_unit_raw
colnames(rgrs_for_splmnt) <- new_colnames


## Apply site aliases
for(i in 1:length(rgrs_for_splmnt$Site)){
  rgrs_for_splmnt$Site[i] <- common_nms[[rgrs_for_splmnt$Site[i]]]
}

## Reorder N to S
rgrs_for_splmnt <- rgrs_for_splmnt[order(match(rgrs_for_splmnt$Site, site_order)),]
rownames(rgrs_for_splmnt) <- NULL
rgrs_for_splmnt <- rgrs_for_splmnt %>% arrange(desc(row_number()))

# Round to 3 significant figures
for(i in new_colnames){
  if(is.numeric(rgrs_for_splmnt[[i]])){
    rgrs_for_splmnt[[i]] <- signif(rgrs_for_splmnt[[i]], digits = 3)}
}

# Burn:Thin -> Thin:Burn; None -> Control
rgrs_for_splmnt[rgrs_for_splmnt$Treatment=='Burn+Thin', 'Treatment'] = 'Thin+Burn'
rgrs_for_splmnt[rgrs_for_splmnt$Treatment == 'None', 'Treatment'] = 'Control'

# Make scatterplots (climate variables) -----------------------------------

# Collapse to one variable per site
climate_by_site <- rgrs_by_unit_raw[, c('Site', 'ppt_by_waterYr', 
                                'cwd_by_waterYr', 'pdsi_by_waterYr')] %>% group_by(Site) %>% 
  summarise(across(everything(), unique))

for(i in 1:length(climate_by_site$Site)){
  climate_by_site$Site[i] <- common_nms[[climate_by_site$Site[i]]]
}


ppt <- make_rgr_pointPlot('ppt_by_waterYr', rgr_tbl = climate_by_site)
cwd_sitelevel <- make_rgr_pointPlot('cwd_by_waterYr', rgr_tbl = climate_by_site)
pdsi <- make_rgr_pointPlot('pdsi_by_waterYr', rgr_tbl = climate_by_site)

print(ppt)
print(cwd_sitelevel)
print(pdsi)

clim <- ggarrange(ppt, cwd_sitelevel, pdsi,
                  ncol = 2, nrow = 2, labels=c('A','B', 'C'), 
                  align = 'hv',
                  common.legend=TRUE, legend = 'bottom')

ggsave(here::here(fig_dir, 'FigS1_clim_boxplots.jpeg'), 
       plot = clim, device = 'jpeg', dpi = 'print',
       width = 6.5, height = 8, units = 'in', scale = 1)


# Make boxplots (rSDI, soils, and intensity) ------------------------------------------

rsdi <- make_rgr_boxplot('rSDI', rgr_tbl = rgrs_for_splmnt)
depth <- make_rgr_boxplot('TotalDepth', rgr_tbl = rgrs_for_splmnt)
awc <- make_rgr_boxplot('AWC_50cm', rgr_tbl = rgrs_for_splmnt)

## For intensity, remove LaTour because pre-treatment measurements are unreliable
intnsty_boxplot <- make_rgr_boxplot('BA_reduction', rgr_tbl = rgrs_for_splmnt[rgrs_for_splmnt$Site != 'LaTour',])
intnsty_boxplot <- intnsty_boxplot + ylab('Treatment intensity \n(% basal area reduction)')

print(rsdi)
print(depth)
print(awc)
print(intnsty_boxplot)

soils <- ggarrange(depth, 
                   awc + ylab('AWC in top 50cm \n(cm/cm)'), 
                   ncol = 2, nrow = 1, labels=c('A','B'), 
                   align = 'hv',
                   common.legend=TRUE, legend = 'bottom')

ggsave(here::here(fig_dir, 'FigS2_soils_boxplots.jpeg'), 
       plot = soils, device = 'jpeg', dpi = 'print',
       width = 6.5, height = 4, units = 'in', scale = 1)

ggsave(here::here(fig_dir, 'FigS3_rsdi_boxplot.jpeg'), 
       plot = rsdi + ylab('Relative stand density \nindex'), device = 'jpeg', dpi = 'print',
       width = 6.5, height = 4, units = 'in', scale = 1)

ggsave(here::here(fig_dir, 'FigS4_intensity_boxplot.jpeg'),
       plot = intnsty_boxplot, device = 'jpeg', dpi = 'print',
       width = 6.5, height = 4, units = 'in', scale = 1)
