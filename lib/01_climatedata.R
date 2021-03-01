
# Packages ----------------------------------------------------------------

library(tidyverse)
library(raster)
library(lubridate)
library(sf)
library(patchwork)

# Studyregion -------------------------------------------------------------

studyregion <- read_sf("data/gis/studyregion_epsg3035.gpkg")

studyregion_latlng <- st_transform(studyregion, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

cntrs <- read_sf("data/admin/countries_europe.shp")

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

dat_sm <- dat_tmp %>%
  as_tibble() %>%
  mutate(gridindex = 1:n()) %>%
  gather(key = key, value = sm, -gridindex) %>%
  filter(!is.na(sm)) %>%
  mutate(year = substring(key, 2, 5),
         month = substring(key, 7, 8)) %>%
  dplyr::select(-key)

# Load ERA5 vapour pressure deficit data ----------------------------------

# Data was downloaded as netCDF from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview

vars <- c("t2m", "d2m")

out <- vector("list", length(vars))

for (k in 1:length(vars)) {
  
  print(vars[k])
  
  dat_ras <- brick("data/climate/era5_temp.nc", var = vars[k])
  dat_ras <- stack(dat_ras)
  
  dat_ras <- crop(dat_ras, studyregion_latlng)
  dat_ras <- raster::mask(dat_ras, studyregion_latlng)
  
  out[[k]] <- dat_ras %>%
    as.matrix()
  
}

# Formula from https://journals.ametsoc.org/view/journals/apme/54/6/jamc-d-14-0321.1.xml?tab_body=pdf
# Based on Allen, R. G., et al. (1998), Crop evapotranspiration. Guidelines for computing crop water requirements. Irrigation and Drainage Paper 56, FAO,Rome, Italy, 27 pp.

vpd_calc <- function(t, td, c1 = 0.611, c2 = 17.67, c3 = 243.5) {
  c1 * exp((c2 * (t - 273.15)) / (t - 273.15 + c3)) - c1 * exp((c2 * (td - 273.15)) / (td - 273.15 + c3))
}

dat_tmp <- out %>% # Calculate VPD
  Reduce(vpd_calc, .) 

rm(out)
gc()

# Bind together into tibble

dat_vpd <- dat_tmp %>%
  as_tibble() %>%
  mutate(gridindex = 1:n()) %>%
  gather(key = key, value = vpd, -gridindex) %>%
  filter(!is.na(vpd)) %>%
  mutate(year = substring(key, 2, 5),
         month = substring(key, 7, 8)) %>%
  dplyr::select(-key)

# Calculate average summer (JJA) SM and VPD --------------------------------

dat <- dat_sm %>%
  left_join(dat_vpd)

write_csv(dat, "data/climate/era5_sm_vpd.csv")
dat <- read_csv("data/climate/era5_sm_vpd.csv")

reference_period <- 1980:2010

dat_summary <- dat %>%
  filter(month %in% c("06", "07", "08")) %>%
  group_by(year, gridindex) %>%
  summarize(sm = mean(sm),
            vpd = mean(vpd)) %>%
  group_by(gridindex) %>%
  mutate(sm_z = (sm - mean(sm[year %in% reference_period])) / sd(sm[year %in% reference_period]),
         vpd_z = (vpd - mean(vpd[year %in% reference_period])) / sd(vpd[year %in% reference_period])) %>%
  ungroup()

write_csv(dat_summary, "data/climate/era5_sm_vpd_summer_anomaly.csv")

# Create reference grid and save as vector --------------------------------

gridindex <- brick("data/climate/era5_sm.nc", var = "swvl1")
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

# Plots -------------------------------------------------------------------

dat_summary_2018 <- dat_summary %>% 
  filter(year %in% c("2018"))

dat_summary_2018_grid <- grid %>% 
  right_join(dat_summary_2018)

st_centroid(dat_summary_2018_grid) %>%
  st_write("data/climate/era5_sm_vpd_2018_gridpoints_epsg3035.gpkg")
  
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, st_crs(cntrs))
st_crs(world) <- st_crs(cntrs)
world <- st_crop(world, st_bbox(cntrs))

p <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = dat_summary_2018_grid, aes(fill = sm_z, col = sm_z)) +
  geom_sf(data = world, color = "black", fill = NA) +
  scale_fill_gradient2(low = "#b2182b", mid = "#FFFFFF", high = "#2166ac",
                       breaks = c(-5, -2.5, 0, 2.5, 5)) +
  scale_color_gradient2(low = "#b2182b", mid = "#FFFFFF", high = "#2166ac",
                        breaks = c(-5, -2.5, 0, 2.5, 5)) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = "Soil moisture anomaly", col = "Soil moisture anomaly") +
  coord_sf(expand = FALSE) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0", color = "black", size = 1.25),
        legend.key.height = unit(0.125, "cm"),
        legend.key.width = unit(0.6, "cm"),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(color = "black", fill = "white"),
        legend.title = element_text(size = 9),
        legend.direction = "horizontal",
        legend.text = element_text(size = 6),
        legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  xlim(st_bbox(world)[1], st_bbox(world)[3]) +
  ylim(st_bbox(world)[2], st_bbox(world)[4]) +
  guides(fill = guide_colorbar(title.position = "top"),
         col = guide_colorbar(title.position = "top"))

ggsave("results/soilmoisture_anomaly.pdf", p, width = 3.5, height = 3.5)

dat_summary_annual <- dat_summary %>%
  group_by(year) %>%
  summarise_at(.vars = c("sm", "vpd"), .funs = mean) %>%
  mutate(sm_z = (sm - mean(sm[year %in% reference_period])) / sd(sm[year %in% reference_period]),
         vpd_z = (vpd - mean(vpd[year %in% reference_period])) / sd(vpd[year %in% reference_period]))
  
p <- ggplot(dat_summary_annual) +
  geom_line(aes(x = as.integer(year), y = sm_z), col = "#BBBBBB") +
  geom_point(aes(x = 2018, y = sm_z[year == 2018])) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0.1, 0.1), "cm")) +
  labs(x = "Year", y = "z-value")



