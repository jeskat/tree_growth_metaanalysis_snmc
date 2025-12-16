###----------------------------------------------------------------------------##
## Create Supplementary Table E1 (Gelman-Rubin statistics for all SSM parameters)
###----------------------------------------------------------------------------##

### Dependencies
library(tidyr)
library(dplyr)
source(here::here('Figures_and_Tables/Code/plotting_config.R'))

## Read GR diagnostics (created by batch_convergence_diag.R)
grdf <- read.csv(file.path(
  here::here('State_space_growth_models/Outputs/Final_SSM_diagnostics/all_gr_diagnostics.csv')))


## Add columns for plotting
grdf$Site <- NULL
grdf$PFT <- NULL
grdf$var_pretty <- NULL

## Get pretty names
for(i in 1:length(grdf$site)){
  grdf$Site[i] <- common_nms[grdf$site[i]]
  grdf$PFT[i] <- pft_dict[grdf$pft[i]]
  grdf$var_pretty[i] <- pretty_nms[grdf$param[i]]
}

# Add a column with point estimate and upper CI
grdf$mean_ci <- paste0(
  round(grdf$Point.est., 2), ' (', round(grdf$Upper.C.I., 2), ')', 
  sep = '')


# Reshape data
gr_table <- spread(grdf[,c('var_pretty', 'mean_ci', 'Site', 'PFT')], 
                   var_pretty, mean_ci)

# Reorder columns
gr_table <-gr_table[,c('Site', 'PFT', ssm_param_order)]

# Order North-South
gr_table <- gr_table[order(match(gr_table$Site, site_order)),]
rownames(gr_table) <- NULL
gr_table <- gr_table %>% arrange(desc(row_number()))


write.csv(apply(gr_table,2,as.character), 
          here::here('Figures_and_Tables/Tables/TabE1_gr_diagnostics.csv'), 
          row.names = FALSE)
