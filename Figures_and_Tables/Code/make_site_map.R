
# Data and dependencies ---------------------------------------------------

library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(ggspatial)
library(dplyr)
library(purrr)

## Load polygons for each unit
plot_poly <- st_read(
  'Covariate_and_metaregressor_processing/Input_data/unit_polygons.shp')

## Load the Sierra Nevada ecoregion boundary
sn_bounds <- st_read(
  'Figures_and_Tables/Map_input_data/Sierra_Nevada_Conservancy_Subregions/Sierra_Nevada_Conservancy_Subregions.shp')

# Dissolve interior polygons
sn_bounds <- st_union(st_buffer(sn_bounds, 0.0001))

## Load elevation raster
elev <- rast('Figures_and_Tables/Map_input_data/SRTM_clippedCA.tif')

## Convert to projected coordinate system (CA Albers)
plots_projected <- st_transform(plot_poly, crs = 3310)

elev_projected <- project(
  x = elev, 
  y = st_crs(plots_projected)$wkt,   
  method = "bilinear"      # Use 'bilinear' or 'cubic' for continuous data like elevation
)

## Find the centroid of each site ------------------------------------------

site_centroids <- plots_projected %>%
  group_by(Site) %>% 
  summarise(geometry = st_union(geometry)) %>% 
  # Use map to calculate the bbox for each individual site geometry
  mutate(bbox = map(geometry, st_bbox)) %>% 
  mutate(
    center_x = map_dbl(bbox, ~ (.x["xmin"] + .x["xmax"]) / 2),
    center_y = map_dbl(bbox, ~ (.x["ymin"] + .x["ymax"]) / 2)
  ) %>%
  # Convert the calculated X/Y into a point geometry
  st_as_sf(coords = c("center_x", "center_y"), crs = st_crs(plots))



# Add experiment attributes -----------------------------------------------

## Make an Experiment columns
site_centroids$Experiment = 'Thin'
site_centroids[site_centroids$Site %in% c('Sequoia', 'TharpsCreek'), 
                'Experiment'] = 'Burn'
site_centroids[site_centroids$Site %in% c('Blodgett', 'STEF', 'Teakettle'), 
                'Experiment'] = 'Thin x Burn'


## Add the number of observations at each site

for(site in c('Blodgett', 'LaTour', 'STEF', 'Teakettle', 'WLakeTahoe')){ #
  tree_data <- read.csv(here::here('State_space_growth_models/Input_data/',
                                   site,
                                   'tree_attrs.csv'))
  if(site %in% c('Blodgett', 'STEF', 'Teakettle')){
    pfts <- tree_data[tree_data$PFT %in% c('cedar', 'fir', 'sugar_pine', 'yellow_pine'),]
  }else if(site %in% c('LaTour', 'WLakeTahoe')){
      pfts <- tree_data[tree_data$PFT %in% c('fir'),]
  }else if(site=='Sequoia'){
    tree_data[tree_data$PFT %in% c('cedar', 'fir', 'sugar_pine'),]
  }else if(site == 'TharpsCreek'){
    tree_data[tree_data$PFT %in% c('fir', 'sugar_pine'),]
  }
  
  site_centroids[site_centroids$Site==site, 'Number of observations'] = sum(pfts$n_obs)
}

p <- ggplot() +
  geom_spatraster(data = elev_projected, alpha=0.5) + 
  scale_fill_wiki_c() + 
  geom_sf(data = sn_bounds, fill = NA) + 
  geom_sf(data = site_centroids, aes(size = `Number of observations`, color = Experiment)) + 
  scale_color_manual(values = c('#D55E00', '#42017C', '#0072B2')) +
  scale_size_continuous(limits = c(0, 30000), 
                        breaks = c(100, 1000, 10000))+
  annotation_scale(location = 'bl') +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         height = unit(1, 'cm'), width = unit(1, 'cm'), 
                         pad_y = unit(1, 'cm'), pad_x = unit(0.5, 'cm')) +
  labs(fill= 'Elevation (m)', 
       color = "Experiment", 
       size = "Number of observations") +
  
  theme(panel.background = element_rect(fill = NA, colour = 'NA')) +
  theme(axis.line = element_line(color='gray'),
        axis.text=element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave('Code/figures/site_map.png', width = 8, height = 6, units ='in')
