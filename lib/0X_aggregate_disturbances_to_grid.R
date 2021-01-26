
# Packages ----------------------------------------------------------------

library(tidyverse)
library(raster)
library(lubridate)
library(sf)
library(patchwork)
library(fasterize)

# Studyregion and grid ----------------------------------------------------

# The following code needs disturbance maps downloaded from 10.5281/zenodo.3924381
path_dist_maps <- "../../Disturbance Europe/Mapping/results/version1.1/"

countries <- list.files(path_dist_maps)

grid <- read_sf("data/climate/climategrid_epsg3035.gpkg")

# Aggregate to climate grid -----------------------------------------------

for (i in 1:length(countries)) {
  
  cntr <- countries[i]
  
  print(cntr)
  
  disturbance <- raster(paste0(path_dist_maps, cntr, "/disturbance_year_", cntr, ".tif"))
  #forest <- raster(paste0(path_dist_maps, cntr, "/forest_cover_", cntr, ".tif"))
  
  ext <- as(extent(disturbance), 'SpatialPolygons')
  ext <- st_as_sf(ext)
  st_crs(ext) <- st_crs(grid)
  
  grid_sel <- st_intersection(st_as_sf(grid), st_as_sf(ext))
  grid_sel_ras <- fasterize(grid_sel, disturbance, field = "gridindex")
  grid_values <- values(grid_sel_ras)
  
  dat <- data.frame(gridindex = grid_values,
                    disturbance = values(disturbance),
                    country = cntr) %>%
    na.omit(.) %>%
    group_by(gridindex, year = disturbance) %>%
    summarize(disturbance_ha = n() * 0.09,
              country = unique(country)) %>%
    ungroup(.)
  
  # forest <- data.frame(gridindex = grid_values,
  #                      forest = values(forest)) %>%
  #   filter(!is.na(forest)) %>%
  #   group_by(gridindex) %>%
  #   summarize(forest_ha = sum(forest == 1, na.rm = TRUE) * 0.09,
  #             land_ha = n() * 0.09) %>%
  #   ungroup(.)
  
  # dat <- dat %>% 
  #   left_join(forest, by = "climate_id")
  
  write_csv(dat, paste0("data/disturbances/aggregated_to_grid/disturbances_aggregated_to_grid_", cntr, ".csv"))
  
}

