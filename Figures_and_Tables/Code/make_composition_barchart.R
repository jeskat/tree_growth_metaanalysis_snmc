###---------------------------------------------------------------------------##
## Makes supplementary figure showing pre-treatment species composition by site
###---------------------------------------------------------------------------##


library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

source(here::here('Figures_and_Tables/Code/plotting_config.R'))

### Load pre-treatment species compositions from each site
## Instantiate list to hold site data
all_sites <- vector(mode = 'list', length = length(site_list))

for(i in 1:length(site_list)){
  site <- site_list[i]
  temp <- read.csv(here::here(
    'Covariate_and_metaregressor_processing/Processed_data/Stand_structure/',
                             paste0(site,'_pretrt_comp.csv')))
  temp$Site <- common_nms[[site]]
  all_sites[[i]] <- temp
}

## Combine into a dataframe
all_sites <- do.call(rbind, all_sites)

## Average species composition across plots in each site
site_totals <- all_sites %>%
  group_by(Site, Species) %>%
  summarise(Species_frac = mean(Species_frac))

## Clean and consolidate species

# Based on the 2020 CFI data, all PIPO/PIJE trees at Latour are PIPO
site_totals[site_totals$Species == 'PIPO/PIJE', 'Species'] = 'PIPO'

# Assign low-incidence species to "other"
site_totals$Species_plot = site_totals$Species
site_totals[!(site_totals$Species %in% c('ABCO', 'CADE', 'PILA', 'QUKE', 'PIJE', 
                                         'ABMA', 'PIPO','PSME', 'SEGI')), 
            'Species_plot'] = 'Other'

# Make species a factor
site_totals$Species_plot = factor(site_totals$Species_plot,
                                  levels = c('ABCO','ABMA','CADE', 'PIJE', 
                                             'PIPO', 'PILA', 'PSME', 'SEGI', 
                                             'QUKE', 'Other'))


pal <- list(
  ABCO = '#88CCEE', 
  ABMA = '#0072B2', 
  CADE = '#117733', 
  PIJE = '#FFEA83', 
  PIPO = '#DDCC77', 
  PILA = '#CC6677', 
  PSME = '#E69F00', 
  SEGI = '#882255', 
  QUKE = '#332288', 
  Other = 'grey')

### Stacked bar plot
barplot <- ggplot(site_totals, aes(fill = Species_plot, y = Species_frac, x = Site)) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = pal) + 
  ylab('Percent basal area (trees >= 15cm DBH)') + 
  labs(fill = 'Species') + 
  theme_minimal(base_size = 14) +
  theme(axis.line = element_line(color='black'), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
  
  
ggsave(here::here(fig_dir, 'FigS4_composition_barchart.jpeg'), 
       plot = barplot, device = 'jpeg', dpi = 'print',
       width = 6.5, height = 5, units = 'in', scale = 1)
