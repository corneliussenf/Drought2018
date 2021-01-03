

# Packages ----------------------------------------------------------------

library(tidyverse)
library(raster)
library(lubridate)
library(sf)

# Studyregion -------------------------------------------------------------

studyregion <- read_sf("data/gis/studyregion_epsg3035.gpkg")

studyregion_latlng <- st_transform(studyregion, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Load ERA5 soil moisture data --------------------------------------------

# Data was downloaded as netCDF from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview

vars <- c("swvl1", "swvl2", "swvl3", "swvl4")

out <- vector("list", length(vars))

for (k in 1:length(vars)) {
  
  print(vars[k])
  
  dat_ras <- brick("data/climate/era5_sm.nc", var = vars[k])
  dat_ras <- stack(dat_ras)
  
  dat_ras <- crop(dat_ras, studyregion_latlng)
  dat_ras <- raster::mask(dat_ras, studyregion_latlng)
  
  out[[k]] <- dat_ras %>%
    as.matrix()
  
}

dat_tmp <- out %>% Reduce("+", .) # Sum soil moisture fomr different layers

rm(out)
gc()

# Bind together into tibble

dat <- dat_tmp %>%
  as_tibble() %>%
  mutate(gridindex = 1:n()) %>%
  gather(key = key, value = sm, -gridindex) %>%
  filter(!is.na(sm)) %>%
  mutate(year = substring(key, 2, 5),
         month = substring(key, 7, 8)) %>%
  dplyr::select(-key)


# Calculate average summer (JJA) soil moisture ----------------------------

dat_summary <- dat %>%
  filter(month %in% c("06", "07", "08")) %>%
  group_by(year, gridindex) %>%
  summarize(sm = mean(sm)) %>%
  group_by(gridindex) %>%
  mutate(sm_z = (sm - mean(sm)) / sd(sm)) %>%
  ungroup()


# Create reference grid and save as vector --------------------------------

gridindex <- brick("data/climate/era5_sm.nc", var = vars[1])
gridindex <- subset(stack(gridindex), 1)
gridindex <- crop(gridindex, studyregion_latlng)
gridindex <- raster::mask(gridindex, studyregion_latlng)
values(gridindex) <- 1:ncell(gridindex)
names(gridindex) <- "gridindex"

grid <- gridindex %>%
  rasterToPolygons() %>% 
  st_as_sf(grid) %>%
  st_transform(., st_crs(studyregion))

grid <- grid %>%
  filter(gridindex %in% dat_summary$gridindex)

write_sf(grid, "data/climate/climategrid_epsg3035.gpkg")

# Add soil moisture to grid, save, and create plot -------------------------

dat_summary_2018 <- dat_summary %>% 
  filter(year %in% c("2018"))

dat_summary_grid <- grid %>% 
  right_join(tmp) %>%
  filter(!is.na(sm_z))

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, st_crs(studyregion))
st_crs(world) <- st_crs(studyregion)
world <- st_crop(world, st_bbox(studyregion))

p <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = dat_summary_grid, aes(fill = sm_z, col = sm_z)) +
  geom_sf(data = world, color = "black", fill = NA) +
  scale_fill_gradient2(low = "#CC3311", mid = "white", high = "#33BBEE", limits = c(-4.5, 4.5)) +
  scale_color_gradient2(low = "#CC3311", mid = "white", high = "#33BBEE", limits = c(-4.5, 4.5)) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = NULL, col = NULL) +
  coord_sf(expand = FALSE) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0"),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.75, "cm"),
        legend.position = "bottom",
        legend.title.align = 0.5,
        legend.direction = "horizontal",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(fill = guide_colorbar(title = "Soil moisture (z-value)\n", 
                               title.position = "bottom",
                               title.align = 0.5),
         color = guide_colorbar(title = "Soil moisture (z-value)\n", 
                               title.position = "bottom",
                               title.align = 0.5))
  # facet_wrap(~year, ncol = 4) +
  # geom_text(data = data.frame(year = 2017:2020), aes(x = Inf, y = Inf, label = year),
  #           hjust = 1.25, vjust = 5, size = 2.5)

ggsave("results/figures/sm_anomaly_2018.pdf", p, width = 2.5, height = 3)
