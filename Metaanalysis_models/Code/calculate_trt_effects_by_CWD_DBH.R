###---------------------------------------------------------------------------##
## Runs meta-analyses of treatment effects on growth for different combinations  
## of PFT, DBH, and CWD. Makes Figures 3, 4, and 5. 
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
library(scico)

source(here::here('Metaanalysis_models/Code/metaanalysis_config.R'))
source(here::here('Metaanalysis_models/Code/metaanalysis_functions.R'))

source(here::here('Figures_and_Tables/Code/plotting_config.R'))
source(here::here('Figures_and_Tables/Code/plotting_functions.R'))

### Load table and list of covariance matrices that contain growth outcomes for 
### all combinations of site, PFT, DBH, and CWD
all_dbh_cwd_results <- read.csv(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/',
                                   'growthOutcomesByCWDAndDBH_', 
                                   date, '.csv'))
all_dbh_cwd_cvMatrices <- readRDS(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/',
                                         'growthOutcomesByCWDAndDBH_', 
                                         date, '.RData'))


###-----------------------------------------------------------------------------
### Run meta-analyses of treatment effects on growth for all combinations of 
### PFT, DBH (15-100cm), and CWD (-3 - 3 standard deviations)
###-----------------------------------------------------------------------------

### Reformat inputs 
## Prep y-values for metafor
dbh_cwd_frmttd <- formatForMetafor(all_dbh_cwd_results)

## Convert list elements to matrices
dbh_cwd_cov_m <- lapply(all_dbh_cwd_cvMatrices, function(x) as.matrix(x))

## Compile treatment effects and p-values
pvals_list <- dbh_cwd_frmttd %>%
  group_by(DBH, CWD) %>% # <<< SPLIT
  summarise(
    # The entire result data frame from test_function is wrapped in list()
    # This creates a "list-column" where each element is a data frame
    results = list({ 
      
      # 1. Calculate the unique key for the current group
      current_key <- paste('DBH', cur_group()$DBH, 'CWD', cur_group()$CWD, sep = '_')
      
      # 2. Look up the corresponding covariance matrices
      # 'dbh_cwd_cov_m' must be available in the environment
      cov_list <- dbh_cwd_cov_m[grepl(current_key, names(dbh_cwd_cov_m))]
      
      # 3. Access the current group's full data (equivalent to .SD)
      group_data <- cur_data()
      
      # 4. Call the test function
      result_df <- get_trt_effects_tbl(full_data = group_data,
                                       cv_list = cov_list,
                                       formula = formula(~Burn*Thin)
      )
      # Return the data frame (it is stored in the list-column)
      result_df
    }),
    .groups = 'drop' # Stop grouping after summarise
  )

## Convert to table
all_dbh_cwd_pvals <- pvals_list %>% 
  tidyr::unnest(cols = results)

### Estimate meta-analysis models for each combination of PFT, DBH, and CWD
## Compile growth outcomes by treatment
all_growth_outcomes <- dbh_cwd_frmttd %>%
  group_by(DBH, CWD) %>% # <<< SPLIT
  summarise(
    # The entire result data frame from test_function is wrapped in list()
    # This creates a "list-column" where each element is a data frame
    results = list({ 
      
      # 1. Calculate the unique key for the current group
      current_key <- paste('DBH', cur_group()$DBH, 'CWD', cur_group()$CWD, sep = '_')
      
      # 2. Look up the corresponding covariance matrices
      # 'dbh_cwd_cov_m' must be available in the environment
      cov_list <- dbh_cwd_cov_m[grepl(current_key, names(dbh_cwd_cov_m))]
      
      # 3. Access the current group's full data 
      group_data <- cur_data()
      
      # 4. Run the meta-analysis
      result_df <- get_growth_predictions(full_data = group_data,
                                          cv_list = cov_list, 
                                          formula = formula(~Burn*Thin)
      )
      # Return the data frame (it is stored in the list-column)
      result_df
    }),
    .groups = 'drop' # Stop grouping after summarise
  )

## Convert to table
all_growth_table <- all_growth_outcomes %>% 
  tidyr::unnest(cols = results)


###-----------------------------------------------------------------------------
### Make a Figure 3 (growth vs. DBH when CWD = 0)
###-----------------------------------------------------------------------------

## Subset results for CWD = 0
pval_tbl_cwd0 <- all_dbh_cwd_pvals[all_dbh_cwd_pvals$CWD==0,]
growth_cwd0 <- all_growth_table[all_growth_table$CWD == 0, ]

## Identify DBH ranges with significant treatment effects
sgnf_dbh <- identify_significant_ranges(pval_tbl_cwd0,DBH)

## Plot and label
dbh_plot <- plot_growth_by_xVar(sgnf_dbh,
                                growth_cwd0, 
                                x_var = DBH)
dbh_plot <- dbh_plot + xlab('DBH (cm)')

print(dbh_plot)

ggsave(here::here('Figures_and_Tables/Figures/Fig3_growth_vs_dbh.png'), 
       plot = dbh_plot, device = 'png',
       width = 6.5, height = 6.5, units = 'in')


###-----------------------------------------------------------------------------
### Make Figure 4 (growth vs. CWD when DBH = 30 cm)
###-----------------------------------------------------------------------------

## Subset results for DBH = 30
pval_tbl_dbh30 <- all_dbh_cwd_pvals[all_dbh_cwd_pvals$DBH==30,]
growth_dbh30 <- all_growth_table[all_growth_table$DBH == 30, ]

## Identify DBH ranges with significant treatment effects
sgnf_cwd <- identify_significant_ranges(pval_tbl_dbh30, CWD)

## Plot
cwd_plot <- plot_growth_by_xVar(sgnf_cwd,
                                growth_dbh30, 
                                x_var = CWD)
cwd_plot <- cwd_plot + xlab(expression("CWD"^T*" (standard deviations from site mean)"))

print(cwd_plot)

ggsave(here::here('Figures_and_Tables/Figures/Fig4_growth_vs_cwd.png'), 
       plot = cwd_plot, device = 'png',
       width = 6.5, height = 6.5, units = 'in')

###-----------------------------------------------------------------------------
### Make Figure 5 (facetted heatmap showing treatment effects for each DBH-CWD 
### combination.
###-----------------------------------------------------------------------------

## Set significance thresholds
all_dbh_cwd_pvals[signif(all_dbh_cwd_pvals$pval, 2) <= 0.05, 'plot_pval'] = 'Strongly significant'
all_dbh_cwd_pvals[signif(all_dbh_cwd_pvals$pval, 2) > 0.05 & 
                    signif(all_dbh_cwd_pvals$pval, 2) <= 0.1, 'plot_pval'] = 'Weakly significant'
all_dbh_cwd_pvals[signif(all_dbh_cwd_pvals$pval, 2) > 0.1, 'plot_pval'] = 'Not significant'

## Convert 'significance' to a factor to control its order
all_dbh_cwd_pvals$plot_pval <- factor(all_dbh_cwd_pvals$plot_pval, levels = c(
  "Strongly significant", "Weakly significant", "Not significant"
))

## Make DBH a character and set its order
all_dbh_cwd_pvals$DBH_char <- factor(as.character(all_dbh_cwd_pvals$DBH), 
                                     levels = as.character(sort(unique(all_dbh_cwd_results$DBH))))

## Make labels legible
all_dbh_cwd_pvals[all_dbh_cwd_pvals$Effect=="Burn:Thin Interaction", "Effect"] = "Burn:Thin\n Interaction"


## Construct heatmap
heatmap <- ggplot(all_dbh_cwd_pvals[all_dbh_cwd_pvals$Effect != "Control",],
                  aes(x = CWD, y = DBH_char, fill = Estimate)) + 
  geom_tile_pattern(
    aes(pattern = plot_pval), # Map significance category to pattern
    pattern_colour  = NA,
    pattern_fill = 'black',
    pattern_density = 0.1,
    # pattern_spacing = 0.025,
    color = "black" # This is the border color of the tiles
    ) +
  scale_fill_scico( # Map color to treatment effect size
    palette = 'bam',
    na.value = 'grey90',
    name = "Mean effect size",
    ## Add limits to make sure color bar is centered on 0
    limits = c(-(max(all_dbh_cwd_pvals$Estimate, na.rm = TRUE)),
               max(all_dbh_cwd_pvals$Estimate, na.rm = TRUE))
    ) +
  # Assign specific patterns to each category
  scale_pattern_manual(
    name = "Significance",
    values = c(
      "Strongly significant" = "crosshatch",
      "Weakly significant" = "stripe",
      "Not significant" = "none" # Use 'none' for no pattern
    )
    ) +
  guides(
    pattern = guide_legend(override.aes = list(
      fill = c("white", "white", "white") 
      ))
    ) +
  labs(
    y = "DBH (cm)",
    x = expression("CWD"^T*" (standard deviations from site mean)")
    ) +
  theme_minimal(base_size = 14) + # A clean theme for the plot
  theme(axis.line = element_line(color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.box = 'vertical') +
  coord_fixed() # Ensures the tiles are square

## Facet by PFT and treatment   
heatmap = heatmap + facet_grid(PFT~factor(Effect, 
                                          levels = c('Burn', 'Thin', 'Burn:Thin\n Interaction')))

print(heatmap)

ggsave(here::here('Figures_and_Tables/Figures/Fig5_faceted_heatmap.png'), 
       plot = heatmap, device = 'png',
       width = 7, height = 8.5, units = 'in')

