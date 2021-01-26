# Packages ----------------------------------------------------------------

library(tidyverse)
library(sf)
library(patchwork)


# Load data ---------------------------------------------------------------

studyregion <- read_sf("data/gis/studyregion_epsg3035.gpkg")

grid <- read_sf("data/climate/climategrid_epsg3035.gpkg")

dist_agg <- list.files("data/disturbances/aggregated_to_grid", pattern = ".csv$", full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows()


# Calculate excess mortality ----------------------------------------------

dat <- dist_agg %>%
  group_by(gridindex, year) %>%
  summarize(disturbance_ha = sum(disturbance_ha, na.rm = TRUE)) %>%
  ungroup() %>%
  split(.$gridindex) %>%
  map(~ right_join(., data.frame(gridindex = unique(.$gridindex),
                                 year = 1986:2020), 
                   by = c("gridindex", "year"))) %>%
  map(~ mutate_at(., .vars = vars(disturbance_ha), .funs = function(x) ifelse(is.na(x), 0, x))) %>%
  bind_rows()

write_csv(dat, "temp/dat.csv")

dat_grid <- grid %>%
  right_join(dat, by = "gridindex")

ggplot() +
  geom_sf(data = dat_grid %>%
            filter(year == 2020), 
          aes(fill = disturbance_ha), col = NA)


