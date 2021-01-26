
# Packages ----------------------------------------------------------------

library(tidyverse)
library(raster)
library(lubridate)
library(sf)
library(patchwork)

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

write_csv(dat, "data/climate/era5_sm_vpn.csv")
dat <- read_csv("data/climate/era5_sm_vpn.csv")

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
world <- st_transform(world, st_crs(studyregion))
st_crs(world) <- st_crs(studyregion)
world <- st_crop(world, st_bbox(studyregion))

p1 <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = dat_summary_2018_grid, aes(fill = sm_z, col = sm_z)) +
  geom_sf(data = world, color = "black", fill = NA) +
  scale_fill_gradient2(low = "#CC3311", mid = "white", high = "#33BBEE", limits = c(-4.5, 4.5)) +
  scale_color_gradient2(low = "#CC3311", mid = "white", high = "#33BBEE", limits = c(-4.5, 4.5)) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = NULL, col = NULL) +
  coord_sf(expand = FALSE) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0"),
        legend.key.height = unit(0.75, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  labs(title = "Summer soil moisture anomaly 2018",
       subtitle = "Reference period 1980-2010")

p2 <- ggplot() +
  geom_sf(data = world, color = NA, fill = "lightgray") +
  geom_sf(data = dat_summary_2018_grid, aes(fill = vpd_z, col = vpd_z)) +
  geom_sf(data = world, color = "black", fill = NA) +
  scale_fill_gradient2(low = "#33BBEE", mid = "white", high = "#CC3311", limits = c(-4.5, 4.5)) +
  scale_color_gradient2(low = "#33BBEE", mid = "white", high = "#CC3311", limits = c(-4.5, 4.5)) +
  theme_linedraw() +
  labs(x = NULL, y = NULL, fill = NULL, col = NULL) +
  coord_sf(expand = FALSE) +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0"),
        legend.key.height = unit(0.75, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  labs(title = "Summer vapor pressure deficit anomaly 2018",
       subtitle = "Reference period 1980-2010")

dat_summary_annual <- dat_summary %>%
  group_by(year) %>%
  summarise_at(.vars = c("sm", "vpd"), .funs = mean) %>%
  mutate(sm_z = (sm - mean(sm[year %in% reference_period])) / sd(sm[year %in% reference_period]),
         vpd_z = (vpd - mean(vpd[year %in% reference_period])) / sd(vpd[year %in% reference_period]))
  
p3 <- ggplot(dat_summary_annual) +
  geom_line(aes(x = as.integer(year), y = sm_z), col = "#BBBBBB") +
  geom_point(aes(x = 2018, y = sm_z[year == 2018])) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0.1, 0.1), "cm")) +
  labs(x = "Year", y = "z-value")

p4 <- ggplot(dat_summary_annual) +
  geom_line(aes(x = as.integer(year), y = vpd_z), col = "#BBBBBB") +
  geom_point(aes(x = 2018, y = vpd_z[year == 2018])) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0.1, 0.1), "cm")) +
  labs(x = "Year", y = "z-value")

# Stitch together

p <- p1 + p2 + p3 + p4 + 
  plot_layout(ncol = 2, 
              widths = c(0.5, 0.5), 
              heights = c(0.85, 0.15))

ggsave("results/figures/climate_figure.pdf", p, width = 7.5, height = 5)

