###---------------------------------------------------------------------------##
## Makes supplementary figure showing distribution of post-treatment PFT sizes 
## at each site. 
###---------------------------------------------------------------------------##


# Data and dependencies ---------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

source(here::here('Figures_and_Tables/Code/plotting_config.R'))
source(here::here('Metaanalysis_models/Code/metaanalysis_functions.R'))


### Get first observation of tree size for each site and format
## Instantiate list to hold site data
all_sizes <- vector(mode = 'list', length = length(site_list))

for(i in 1:length(site_list)){
  site <- site_list[i]
  ## Get sizes
  tree_obs <- read.csv(here::here(
    'State_space_growth_models/Input_data/', site,
    'dbh_tree_obs.csv'))
  ## Get PFTs
  tree_attrs <- read.csv(here::here(
    'State_space_growth_models/Input_data/', site,
    'tree_attrs.csv'))
  ## Merge
  temp <- inner_join(tree_obs[,1:2], tree_attrs, by = 'TreeID')
  colnames(temp)[2] <- 'DBH'
  temp$Site <- common_nms[[site]]
  all_sizes[[i]] <- temp
}

## Combine into a dataframe
all_sizes <- do.call(rbind, all_sizes)

pft_sizes <- all_sizes[all_sizes$PFT %in% c('cedar', 'fir', 
                                            'white_pine', 'yellow_pine'),]

## Change PFT names
for(p in unique(pft_sizes$PFT)){
  pft_sizes[pft_sizes$PFT==p, 'PFT'] = pft_dict[[p]]
}


# Table of PFT counts -----------------------------------------------------

pft_counts <- pft_sizes %>%
  group_by(Site, PFT) %>%
  summarise(`Number of trees` = n(),
            `Number of growth observations` = sum(n_obs-1)) 

## Order alphabetically
pft_counts <- pft_counts[order(match(pft_counts$Site, site_list)),]


write_csv(pft_counts, here::here(tbl_dir, 'TabS1_pftCounts.csv'))


# PFT size distributions --------------------------------------------------

## Drop site-PFT combinations that we didn't model
sizes_filtered <- pft_sizes[(pft_sizes$Site %in% c('Blodgett', 'STEF', 'Teakettle')) |
                          (pft_sizes$Site %in% c('LaTour', 'W. Lake Tahoe') & 
                             pft_sizes$PFT=='Fir') |
                          (pft_sizes$Site=='Sequoia' & 
                             pft_sizes$PFT %in% c('Fir', 'Cedar', 'Sugar Pine')) |
                          (pft_sizes$Site == "Tharp's Creek" & 
                             pft_sizes$PFT %in% c('Fir', 'Sugar Pine')),]

## Order north to south
sizes_filtered$Site <- factor(sizes_filtered$Site, site_order)

## Grouped boxplots
size_dist <- ggplot(sizes_filtered, 
       aes(x=Site, y=DBH, fill=PFT)) + 
  geom_boxplot() + 
  scale_fill_manual(values=unlist(pft_cols), name='PFT') +
  coord_flip() + theme_minimal(base_size = 14) +
  xlab('Site') + 
  ylab('DBH (cm)') + 
  theme(axis.line = element_line(color='black'), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'bottom')

print(size_dist)

ggsave(here::here(fig_dir, 'FigS5_pft_size_distribution.jpeg'), 
       plot = size_dist, device = 'jpeg', dpi = 'print',
       width = 6.5, height = 5, units = 'in', scale = 1)

## Print additional statistics
median(sizes_filtered$DBH)
