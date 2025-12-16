###----------------------------------------------------------------------------##
## Create a table summarizing parameter estimates from each growth model
###----------------------------------------------------------------------------##


library(tidyr)
library(ggplot2)
library(ggpubr)
library(latex2exp)

source(here::here('Figures_and_Tables/Code/plotting_config.R'))
source(here::here('Figures_and_Tables/Code/plotting_functions.R'))
source(here::here('Metaanalysis_models/Code/metaanalysis_functions.R'))


## All parameters from state-space growth models 
allParams <- read.csv(file.path(here::here('State_space_growth_models/Outputs',  
                                           paste0('allParams_', 
                                                  '2025-12-02', ## Replace date if needed
                                                  '.csv'))))
## Add formatted names
allParams <- addPrettyNames(allParams)


###----------------------------------------------------------------------------##
## Create a Supplementary Tables E2-E3 summarizing parameter estimates from 
## each growth model
###----------------------------------------------------------------------------##

## Remove unnecessary columns
ssm_outs <- allParams[allParams$variable %in% c(fxdEff, othrPrms),]

## Add a column with mean and confidence interval
ssm_outs$mean_ci <- paste0(
  signif(ssm_outs$mean, 2), '\n[', 
  signif(ssm_outs$low_ci, 2), ', ', 
  signif(ssm_outs$high_ci, 2), ']', 
  sep = '')

## Reshape data
ssm_outs$var_pretty <- NaN
for(i in 1:length(ssm_outs$site)){
  ssm_outs$var_pretty[i] <- pretty_nms[ssm_outs$variable[i]]
}

mean_ci_tbl <- spread(ssm_outs[,c('var_pretty', 'mean_ci', 'common_nm', 'PFT')], 
                      var_pretty, 
                      mean_ci)

## Reverse the order (should read North-South)
mean_ci_tbl <- mean_ci_tbl[order(match(mean_ci_tbl$common_nm, site_order), 
                                 decreasing = TRUE),]


## Reorder columns
mean_ci_ssmParams <-mean_ci_tbl[,c('common_nm', 'PFT', ssm_param_order)]

# Save
write.csv(apply(mean_ci_ssmParams,2,as.character), 
          here::here('Figures_and_Tables/Tables/TabE2-E3_mean_ci_ssmParams.csv'),
          row.names = FALSE)


###----------------------------------------------------------------------------##
## Create Supplementary Figure E1 summarizing parameter estimates from each 
## growth model
###----------------------------------------------------------------------------##

### Make each panel

intOverall <- plotSSMparam(allParams, 'int_overall', plot_metanalysis = FALSE)
print(intOverall)

## Harmonize names of the log(DBH) parameter (different across SSMs)
allParams[allParams$variable %in% c('slope_size', 'log_size'), 'variable'] = 'ln_size'

lnSize <- plotSSMparam(allParams, 'ln_size', plot_metanalysis = FALSE)
print(lnSize)

sqSize <- plotSSMparam(allParams, 'd2_size')
print(sqSize)

cwd <- plotSSMparam(allParams, 'betaCWD')
print(cwd)

burn <- plotSSMparam(allParams, 'treat_mean[Burn]', plot_metanalysis = FALSE)
print(burn)

thin <- plotSSMparam(allParams, 'treat_mean[Thin]', plot_metanalysis = FALSE)
print(thin)

burnThin <- plotSSMparam(allParams, 'treat_mean[Burn+Thin]', plot_metanalysis = FALSE)
print(burnThin)

thinCWD <- plotSSMparam(allParams, 'treat_cwd[Thin]')
print(thinCWD)

burnCWD <- plotSSMparam(allParams, 'treat_cwd[Burn]')
print(burnCWD)

thinBurnCWD <- plotSSMparam(allParams, 'treat_cwd[Burn+Thin]')
print(thinBurnCWD)

thinSize <- plotSSMparam(allParams, 'treat_size[Thin]')
print(thinSize)

burnSize <- plotSSMparam(allParams, 'treat_size[Burn]')
print(burnSize)

thinBurnSize <- plotSSMparam(allParams, 'treat_size[Burn+Thin]')
print(thinBurnSize)


plotSD <- plotSSMparam(allParams, 'int_plot_sd', plot_metanalysis = FALSE)
print(plotSD)

treeSD <- plotSSMparam(allParams, 'int_tree_sd', plot_metanalysis = FALSE)
print(treeSD)

obsSD <- plotSSMparam(allParams, 'obs_sd', plot_metanalysis = FALSE)
print(obsSD)

resSD <- plotSSMparam(allParams, 'res_sd', plot_metanalysis = FALSE)
print(resSD)

## Combine into a summary graphic
## First six panels
a_f <- ggarrange(intOverall, lnSize, sqSize, cwd, burn, thin,   
          ncol=2, nrow=3, common.legend=TRUE, legend = 'bottom', align = "hv",
          labels = c('A', 'B', 'C', 'D', 'E', 'F'))

## Second six panels
g_l <- ggarrange(burnThin, thinCWD, burnCWD, thinBurnCWD, thinSize, burnSize, 
                 ncol=2, nrow=3, common.legend=TRUE, legend = 'bottom', align = "hv",
                 labels = c('G', 'H', 'I', 'J', 'K', 'L'))

## Last five panels
m_q <- ggarrange(thinBurnSize, treeSD, plotSD, obsSD, resSD,
                 ncol=2, nrow=3, common.legend=TRUE, legend = 'bottom', align = "hv",
                 labels = c('M', 'N', 'O', 'P', 'Q'))

ggsave(here::here('Figures_and_Tables/Figures/FigE1_SSM_params_A-F.png'),
       plot = a_f, device = 'png',
       width = 7, height = 8.5, units = 'in')

ggsave(here::here('Figures_and_Tables/Figures/FigE1_SSM_params_G-L.png'),
       plot = g_l, device = 'png',
       width = 7, height = 8.5, units = 'in')

ggsave(here::here('Figures_and_Tables/Figures/FigE1_SSM_params_M_Q.png'),
       plot = m_q, device = 'png',
       width = 7, height = 8.5, units = 'in')
