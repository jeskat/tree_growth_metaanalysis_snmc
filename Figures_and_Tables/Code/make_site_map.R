
# Data and dependencies ---------------------------------------------------

library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(ggspatial)
library(dplyr)
library(purrr)

source(here::here('Figures_and_Tables/Code/plotting_config.R'))

## Load polygons for each unit
plot_poly <- st_read(
  'Covariate_and_metaregressor_processing/Input_data/unit_polygons.shp')

### Load the Sierra Nevada ecoregion boundary
## Source: Sierra Nevada Conservancy. 2024. “Sierra Nevada Conservancy Subregions.” 
## https://gis.data.cnra.ca.gov/datasets/SNC::sierra-nevada-conservancy-subregions/about.

sn_bounds <- st_read(
  'Figures_and_Tables/Map_input_data/Sierra_Nevada_Conservancy_Subregions/Sierra_Nevada_Conservancy_Subregions.shp')

# Dissolve interior polygons
sn_bounds <- st_union(st_buffer(sn_bounds, 0.0001))

### Load elevation raster
## Source: Jarvis, A., H.I. Reuter, A. Nelson, and E. Guevara. 2008. 
## Hole-Filled Seamless SRTM Data V4. International Centre for Tropical 
## Agriculture (CIAT). https://srtm.csi.cgiar.org.
## The original data has been clipped to the boundaries of CA. 

elev <- rast('Figures_and_Tables/Map_input_data/SRTM_clippedCA.tif')

## Convert to projected coordinate system (CA Albers)
plots_projected <- st_transform(plot_poly, crs = 3310)

elev_projected <- project(
  x = elev, 
  y = st_crs(plots_projected)$wkt,   
  method = "bilinear"      # Use 'bilinear' or 'cubic' for continuous data like elevation
)

# Find the centroid of each site ------------------------------------------

site_centroids <- plots_projected %>%
  group_by(Site) %>% 
  summarise(geometry = st_union(geometry)) %>% 
  # Use map to calculate the bbox for each individual site geometry
  mutate(bbox = map(geometry, st_bbox)) %>% 
  mutate(
    center_x = map_dbl(bbox, ~ (.x["xmin"] + .x["xmax"]) / 2),
    center_y = map_dbl(bbox, ~ (.x["ymin"] + .x["ymax"]) / 2)
  ) %>%
  # Remove the MultiPolygon geometry column
  st_drop_geometry() %>% 
  # Remove the bbox list column (optional, for a cleaner table)
  select(-bbox) %>% 
  # Convert the calculated X/Y into a point geometry
  st_as_sf(coords = c("center_x", "center_y"), crs = st_crs(plots_projected))


# Add experiment attributes -----------------------------------------------

## Make an Experiment columns
site_centroids$Experiment = 'Thin'
site_centroids[site_centroids$Site %in% c('Sequoia', 'TharpsCreek'), 
                'Experiment'] = 'Burn'
site_centroids[site_centroids$Site %in% c('Blodgett', 'STEF', 'Teakettle'), 
                'Experiment'] = 'Thin x Burn'


## Add the number of observations at each site

for(site in site_centroids$Site){ #
  tree_data <- read.csv(here::here('State_space_growth_models/Input_data/',
                                   site,
                                   'tree_attrs.csv'))
  if(site %in% c('Blodgett', 'STEF', 'Teakettle')){
    pfts <- tree_data[tree_data$PFT %in% c('cedar', 'fir', 'sugar_pine', 'yellow_pine'),]
  }else if(site %in% c('LaTour', 'WLakeTahoe')){
      pfts <- tree_data[tree_data$PFT %in% c('fir'),]
  }else if(site=='Sequoia'){
    pfts <- tree_data[tree_data$PFT %in% c('cedar', 'fir', 'sugar_pine'),]
  }else if(site == 'TharpsCreek'){
    pfts <- tree_data[tree_data$PFT %in% c('fir', 'sugar_pine'),]
  }
  
  site_centroids[site_centroids$Site==site, 'Number of observations'] = sum(pfts$n_obs)
}


# Make map ----------------------------------------------------------------

p <- ggplot() +
  geom_spatraster(data = elev_projected, alpha=0.5) + 
  scale_fill_wiki_c() + 
  geom_sf(data = sn_bounds, fill = NA) + 
  geom_sf(data = site_centroids, 
          aes(size = `Number of observations`, color = Experiment)) + 
  scale_color_manual(values = unlist(trt_cols)) +
  scale_size_continuous(limits = c(0, 37000), 
                        breaks = c(100, 1000, 10000))+
  annotation_scale(location = 'bl') +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         height = unit(1, 'cm'), width = unit(1, 'cm'), 
                         pad_y = unit(1, 'cm'), pad_x = unit(0.5, 'cm')) +
  labs(color = "Experiment", 
       size = "Number of \nobservations",
       fill= 'Elevation (m)') +
  
  guides(color = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  
  theme(panel.background = element_rect(fill = NA, colour = 'NA')) +
  theme(axis.line = element_line(color='gray'),
        axis.text=element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(here::here(fig_dir, 'Fig1_map.jpeg'), 
       device = 'jpeg', dpi = 'print',
       width = 6.5, height = 4.88, units ='in')
