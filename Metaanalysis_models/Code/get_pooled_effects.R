###---------------------------------------------------------------------------##
## Runs meta-analyses and makes Figure 2 and Supplementary Tables F1-F3 
## (Pooled effects of treatment, tree size, measurement period climatic water 
## deficit (CWDT), and their interactions on tree growth.)
###---------------------------------------------------------------------------##

###-----------------------------------------------------------------------------
### Dependencies
###-----------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggpattern)
library(latex2exp)
library(tidyverse)
library(metafor)
library(Matrix)
library(clubSandwich)

source(here::here('Metaanalysis_models/Code/metaanalysis_config.R'))
source(here::here('Metaanalysis_models/Code/metaanalysis_functions.R'))

source(here::here('Figures_and_Tables/Code/plotting_config.R'))
source(here::here('Figures_and_Tables/Code/plotting_functions.R'))

## Load data
## Growth effects by unit (reparameterized to DBH = 30cm, CWD = 0)
growth_effects <- read.csv(here::here('Metaanalysis_models/Reparameterized_SSM_outputs/', 
                                      paste0('unit_growth_effects_by_trt_',
                                             date, 
                                             '.csv')))

## Covariance matrices associated with fixed_effects (unit-level)
growth_cov_unit <- readRDS(here::here('Metaanalysis_models/Reparameterized_SSM_outputs/',
                                      paste0('cov_unit_growth_by_trt_', 
                                             date, 
                                             '.RData')))


## Reparameterized growth, log(DBH), and CWD outcomes by treatment (site-level)
fixed_effects <- read.csv(file.path(here::here('Metaanalysis_models/Reparameterized_SSM_outputs/', 
                                               paste0('fixed_effects_by_trt_', 
                                                      date, 
                                                      '.csv'))))

## Covariance matrices associated with fixed_effects (site-level)
size_cov <- readRDS(here::here('Metaanalysis_models/Reparameterized_SSM_outputs/',
                               paste0('cov_size_by_trt_', date, '.RData')))

cwd_cov_by_trt <- readRDS(here::here('Metaanalysis_models/Reparameterized_SSM_outputs/',
                                     paste0('cov_cwd_by_trt_', date, '.RData')))


###-----------------------------------------------------------------------------
### Meta-analyses of treatment, tree size, and measurement period CWD effects 
### on growth
###-----------------------------------------------------------------------------

### Meta-analysis of treatment effects for unit-level data (30cm tree)
## Format inputs
growth_effects$Treatment <- sub(".*\\[([^][]+)].*", "\\1", 
                                growth_effects$variable)
unit_effects <- formatForMetafor(growth_effects)
unit_cov_m <- lapply(growth_cov_unit, function(x) as.matrix(x))

## Run meta-analysis for all PFTs
meta_unit <- make_meta_df(full_data = unit_effects, cv_list = unit_cov_m,
                          formula = formula(~Burn*Thin), by_unit = TRUE)

### Meta-analysis of log(DBH) effect across treatments
## Format inputs
fixed_effects$Treatment <- sub(".*\\[([^][]+)].*", "\\1", 
                               fixed_effects$variable)

size_effects <- formatForMetafor(
  fixed_effects[grepl('size', fixed_effects$variable),]
)
size_cov_m <- lapply(size_cov, function(x) as.matrix(x))
meta_size <- make_meta_df(full_data = size_effects, cv_list = size_cov_m,
                          formula = formula(~Burn*Thin))

## Meta-analysis of CWD effect across treatments
cwd_effects <- formatForMetafor(
  fixed_effects[grepl('cwd', fixed_effects$variable),]
)
cwd_cov_by_trt_m <- lapply(cwd_cov_by_trt, function(x) as.matrix(x))
meta_cwd <- make_meta_df(full_data = cwd_effects, cv_list = cwd_cov_by_trt_m,
                         formula = formula(~Burn*Thin))

###-----------------------------------------------------------------------------
### Construct summary figure (Figure 2)
###-----------------------------------------------------------------------------

## Panel A
cntrlGrowth <- plot_meta_analysis(meta_unit[[1]], 
                                  axis_label = 'Growth (log DBH increment)',
                                  control = TRUE,
                                  zero_line = FALSE)

## Panel B
trtGrowth <- plot_meta_analysis(meta_unit[[1]], 
                                axis_label = 'Effect of treatment on growth',
                                control = FALSE,
                                zero_line = TRUE)

## Combine Panels A & B
growthMeta <- ggarrange(cntrlGrowth, trtGrowth, 
                        nrow = 1,
                        common.legend = TRUE, legend = 'none',
                        labels = c('A', 'B'),
                        widths = c(0.55, 1))

## Panel C
cntrlSize <- plot_meta_analysis(meta_size[[1]], 
                                axis_label = 'Effect of log DBH on growth',
                                control = TRUE,
                                zero_line = TRUE)

## Panel D
trtSize <- plot_meta_analysis(meta_size[[1]], 
                              axis_label = 'Effect of treatment on the log DBH effect',
                              control = FALSE,
                              zero_line = TRUE)

## Combine Panels C & D
sizeMeta <- ggarrange(cntrlSize, trtSize, 
                      nrow = 1,
                      common.legend = TRUE, legend = 'none',
                      labels = c('C', 'D'),
                      widths = c(0.55, 1))

## Panel E
cntrlCWD <- plot_meta_analysis(meta_cwd[[1]], 
                               axis_label = expression("Effect of CWD"^T*" on growth"),
                               control = TRUE,
                               zero_line = TRUE)

## Panel F
trtCWD <- plot_meta_analysis(meta_cwd[[1]], 
                             axis_label = expression("Effect of treatment on the CWD"^T*" effect"),
                             control = FALSE,
                             zero_line = TRUE)

## Combine Panels E & F
cwdMeta <- ggarrange(cntrlCWD, 
                     trtCWD, #+ scale_y_continuous(breaks=seq(-0.3,0.3,0.3)), 
                     nrow = 1,
                     common.legend = TRUE, legend = 'none',
                     labels = c('E', 'F'),
                     widths = c(0.55, 1))

## Combine all three rows
fig2 <- ggarrange(growthMeta, sizeMeta, cwdMeta, 
                  nrow = 3, 
                  common.legend = TRUE, legend = 'bottom', 
                  legend.grob = get_legend(cntrlGrowth))

print(fig2)

ggsave(here::here('Figures_and_Tables/Figures/Fig2_pooled_effects.png'), plot = fig2, device = 'png',
       width = 8, height = 8, units = 'in')

###-----------------------------------------------------------------------------
### Construct summary tables (Supplementary Tables F1-F3)
###-----------------------------------------------------------------------------

## Supplementary Table F1: pooled effects of treatment on growth
pooled_effects_table <- make_results_table(meta_unit[[2]], by_unit = TRUE)

## Supplementary Table F2: pooled effects of treatment on log(DBH)
size_effects_table <- make_results_table(meta_size[[2]], by_unit = FALSE)

## Supplementary Table F3: pooled effects of treatment on measurement period CWD
cwd_effects_table <- make_results_table(meta_cwd[[2]], by_unit = FALSE)

# Save
write.csv(pooled_effects_table, 
          here::here('Figures_and_Tables/Tables/TabF1_pooled_effects_table.csv'),
          row.names = FALSE)
write.csv(size_effects_table, 
          here::here('Figures_and_Tables/Tables/TabF2_pooled_size_table.csv'), 
          row.names = FALSE)
write.csv(cwd_effects_table, 
          here::here('Figures_and_Tables/Tables/TabF3_pooled_cwd_table.csv'),
          row.names = FALSE)
