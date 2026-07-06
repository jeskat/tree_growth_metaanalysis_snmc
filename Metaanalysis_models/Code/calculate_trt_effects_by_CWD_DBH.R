###---------------------------------------------------------------------------##
## Runs meta-analyses of treatment effects on growth for different combinations  
## of PFT, DBH, and CWD. Makes Figures 4, 5, S10, and S11. 
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
                                   'growthOutcomesByCWDAndDBH.csv'))
all_dbh_cwd_cvMatrices <- readRDS(paste0('Metaanalysis_models/Reparameterized_SSM_outputs/',
                                         'growthOutcomesByCWDAndDBH.RData'))


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
### Make Figure S9 (growth vs. DBH when CWD = 0)
###-----------------------------------------------------------------------------

## Subset results for CWD = 0
growth_cwd0 <- all_growth_table[all_growth_table$CWD == 0, ]

## Define custom limits for y-axis for each PFT
dbh_lims <- data.frame(
  PFT = c('Cedar', 'Cedar', 'Fir', 'Fir', 'Sugar Pine', 'Sugar Pine'
          , 'Yellow Pine', 'Yellow Pine'),
  pred = c(0, 0.6,
           0, 0.6,
           0, 1.1
           ,0, 3.5
           )
)

dbh_splmnt <- plot_growth_by_xVar(growth_cwd0, 
                                  DBH, x_subset = c(30, 60, 90), 
                                  logscale = TRUE, dodge_width = 10,
                                  dummy_limits = dbh_lims)

dbh_splmnt <- dbh_splmnt + xlab('DBH (cm) \n')

ggsave(here::here(fig_dir, 'FigS10_growth_vs_dbh.tif'), 
       plot = dbh_splmnt, device = 'tiff', dpi = 600,
       width = 16.5, height = 16.5, units = 'cm')


###-----------------------------------------------------------------------------
### Make Figure S10 (growth vs. CWD when DBH = 30 cm)
###-----------------------------------------------------------------------------

## Subset results for DBH = 30
growth_dbh30 <- all_growth_table[all_growth_table$DBH == 30, ]


cwd_lims <- data.frame(
  PFT = c('Cedar', 'Cedar', 'Fir', 'Fir', 'Sugar Pine', 'Sugar Pine'
          , 'Yellow Pine', 'Yellow Pine'),
  pred = c(0, 2,
           0, 0.8,
           0, 1.3
           ,0, 3.5
           )
)

cwd_splmnt <- plot_growth_by_xVar(growth_dbh30, 
                                CWD, x_subset = c(-3, -2, -1, 0, 1, 2, 3), 
                                logscale = TRUE, dodge_width = 0.5, 
                                dummy_limits = cwd_lims)

cwd_splmnt <- cwd_splmnt + xlab("Measurement period CWD \n(standard deviations from site mean)")

ggsave(here::here(fig_dir, 'FigS11_growth_vs_cwd.tif'), 
       plot = cwd_splmnt, device = 'tiff', dpi = 600,
       width = 16.5, height = 16.5, units = 'cm')


###-----------------------------------------------------------------------------
### Make Figure 4 (unlogged growth vs. DBH and CWD for cedar, fir, and sugar pine)
###-----------------------------------------------------------------------------

## DBH portion
dbh_plot <- plot_growth_by_xVar(growth_cwd0[growth_cwd0$PFT != 'Yellow Pine',], 
                                DBH, x_subset = c(30, 60, 90), 
                                logscale = FALSE, dodge_width = 10, ncol = 1,
                                dummy_limits = dbh_lims[dbh_lims$PFT != 'Yellow Pine',])

dbh_plot <- dbh_plot + xlab('DBH (cm) \n')

## CWD portion
cwd_plot <- plot_growth_by_xVar(growth_dbh30[growth_dbh30$PFT != 'Yellow Pine',], 
                                CWD, x_subset = c(-3, -2, -1, 0, 1, 2, 3), 
                                logscale = FALSE, dodge_width = 0.5, ncol = 1,
                                dummy_limits = cwd_lims[cwd_lims$PFT != 'Yellow Pine',])

cwd_plot <- cwd_plot + xlab("pCWD (standard deviations \nfrom site mean)")

## Combine panels
dbh_and_cwd <- ggarrange(dbh_plot, 
          cwd_plot + theme(axis.title.y = element_blank()),
          nrow = 1,
          common.legend = TRUE,
          legend = 'bottom',
          labels = c('A', 'B'))

dbh_and_cwd <- annotate_figure(
  dbh_and_cwd,
  top = text_grob('Figure 4', size = 12)
)

ggsave(here::here(fig_dir, 'Fig4_growth_vs_dbh_cwd.tif'), 
       plot = dbh_and_cwd, device = 'tiff', dpi = 600,
       width = 18, height = 22, units = 'cm')


###-----------------------------------------------------------------------------
### Make Figure 5 (facetted heatmap showing treatment effects for each DBH-CWD 
### combination.
###-----------------------------------------------------------------------------

## Set flag for using log or unlogged scale
log_scale = TRUE

## Set significance thresholds
all_dbh_cwd_pvals[signif(all_dbh_cwd_pvals$pval, 2) <= 0.05, 'plot_pval'] = "Clear (p <= 0.05)"
all_dbh_cwd_pvals[signif(all_dbh_cwd_pvals$pval, 2) > 0.05 & 
                    signif(all_dbh_cwd_pvals$pval, 2) <= 0.1, 'plot_pval'] = 'Suggestive (0.05 < p <= 0.1)'
all_dbh_cwd_pvals[signif(all_dbh_cwd_pvals$pval, 2) > 0.1, 'plot_pval'] = 'Unclear (p > 0.1)'

## Convert 'significance' to a factor to control its order
all_dbh_cwd_pvals$plot_pval <- factor(all_dbh_cwd_pvals$plot_pval, levels = c(
  "Clear (p <= 0.05)", "Suggestive (0.05 < p <= 0.1)", "Unclear (p > 0.1)"
))

## Make DBH a character and set its order
all_dbh_cwd_pvals$DBH_char <- factor(as.character(all_dbh_cwd_pvals$DBH), 
                                     levels = as.character(sort(unique(all_dbh_cwd_results$DBH))))

## Make labels legible
all_dbh_cwd_pvals[all_dbh_cwd_pvals$Effect=="Thin:Burn Interaction", "Effect"] = "Thin:Burn\n Interaction"

to_plot <- all_dbh_cwd_pvals


## Construct heatmap
heatmap <- ggplot(to_plot[to_plot$Effect != "Control",], 
                  aes(x = CWD, y = DBH_char, fill = Estimate)) + 
  geom_tile_pattern(
    aes(pattern = plot_pval), # Map significance category to pattern
    pattern_colour  = NA,
    pattern_fill = 'black',
    pattern_density = 0.1,
    # pattern_spacing = 0.025,
    color = "black" # This is the border color of the tiles
    ) +
  
  # Assign specific patterns to each category
  scale_pattern_manual(
    name = "Statistical clarity",
    values = c(
      "Clear (p <= 0.05)" = "crosshatch",
      "Suggestive (0.05 < p <= 0.1)" = "stripe",
      "Unclear (p > 0.1)" = "none" # Use 'none' for no pattern
    )
  ) +
  guides(
    pattern = guide_legend(override.aes = list(
      fill = c("white", "white", "white") 
    ))
  ) +
  labs(
    y = "DBH (cm)",
    x = "pCWD (standard deviations \nfrom site mean)"
  ) +
  theme_minimal(base_size = 14) + # A clean theme for the plot
  theme(axis.line = element_line(color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.box = 'vertical') +
  coord_fixed() # Ensures the tiles are square
  
if(log_scale==FALSE){
  heatmap <- heatmap + scale_fill_scico( # Map color to treatment effect size
    palette = 'bam',
    na.value = 'grey90',
    name = "Mean effect size",
    ## Add limits to make sure color bar is centered on 0
    limits = c(-(max(to_plot$Estimate, na.rm = TRUE)),
               max(to_plot$Estimate, na.rm = TRUE)),
    breaks = c(log(0.25), log(0.5), log(1), log(2), log(4)),, # Tells ggplot where to put the tick marks
    labels = c("0.25", "0.5", "1", "2", "4") # Displays the original unlogged numbers
  ) 
  }else{
    heatmap <- heatmap + scale_fill_scico( # Map color to treatment effect size
        palette = 'bam',
        na.value = 'grey90',
        name = "Mean effect size",
        ## Add limits to make sure color bar is centered on 0
        limits = c(-(max(to_plot$Estimate, na.rm = TRUE)),
                   max(to_plot$Estimate, na.rm = TRUE))
      )
        }
  

## Facet by PFT and treatment   
heatmap = heatmap + facet_grid(PFT~factor(Effect, 
                                          levels = c('Burn', 'Thin', 'Thin:Burn\n Interaction')))

heatmap <- annotate_figure(
  heatmap,
  top = text_grob('Figure 5', size = 12)
)
print(heatmap)

ggsave(here::here(fig_dir, 'Fig5_faceted_heatmap.tiff'), 
       plot = heatmap, device = 'tif', dpi = 600,
       width = 18, height = 22, units = 'cm')

