###---------------------------------------------------------------------------##
## Makes supplementary figure showing correlations among metaregressors
###---------------------------------------------------------------------------##


library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(corrplot)

source(here::here('Figures_and_Tables/Code/plotting_config.R'))
source(here::here('Metaanalysis_models/Code/metaanalysis_config.R'))

### Load meta-regressors
rgrs_by_unit_raw <- read.csv(
  here::here('Covariate_and_metaregressor_processing/Processed_data/all_regressors_by_unit.csv'))

### Remove columns we don't use
rgrs <- rgrs_by_unit_raw[,colnames(rgrs_by_unit_raw) %in% 
                           c('Site', 'UnitID', metaregressors)]

## Apply shorthand for names of variables in correlation plot
corrplot_colnames <- colnames(rgrs)
for(i in 1:length(corrplot_colnames)){
  if(colnames(rgrs[i]) %in% names(corrplot_shorthand))
  corrplot_colnames[i] <- corrplot_shorthand[[colnames(rgrs[i])]]
}

colnames(rgrs) <- corrplot_colnames

## Calculate correlations
corr_all <- cor(rgrs[,3:(length(rgrs))])

png(here::here(fig_dir, 'FigS8_metaregressor_correlation.jpeg'), 
    width = 8000 , height = 8000,
    units = "px", type = "quartz", pointsize = 20, res = 800)
corrplot(corr_all, method = 'color', type='upper', addCoef.col = 'black', 
         diag = FALSE)
dev.off()
