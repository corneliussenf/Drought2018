
# Packages ----------------------------------------------------------------

library(tidyverse)
library(sf)
library(patchwork)

# Load data ---------------------------------------------------------------

grid <- read_sf("data/climate/climategrid_epsg3035.gpkg")

cntrs <- read_sf("data/admin/countries_europe.shp")

#grid_eur <- st_join(grid, cntrs, join = st_intersects)
#write_sf(grid_eur, "temp/grid_eur.gpkg")
grid_eur <- read_sf("temp/grid_eur.gpkg")

cntr_groups <- read_csv("data/admin/country_grouping.csv")

dist_agg <- list.files("data/disturbances/aggregated_to_grid", 
                       glob2rx("*disturbance*.csv$"), full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows()

# Drop the year 2020 fpr Norway, as there are many false positives due to data inconsistencies
dist_agg <- dist_agg %>%
  filter(!(country == "norway" & year == 2020))

forest_agg <- list.files("data/disturbances/aggregated_to_grid", 
                         pattern = glob2rx("*forest*.csv$"), full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows()

climate_2018 <- read_csv("data/climate/era5_sm_vpd_summer_anomaly.csv")

# Derive totals -----------------------------------------------------------

dist_agg_total <- dist_agg %>%
  group_by(year) %>%
  summarize(disturbance_ha = sum(disturbance_ha, na.rm = TRUE))

write_csv(dist_agg_total, "results/disturbance_area_map_based_total_europe.csv")

# By drought area

ggplot() +
  geom_line(data = dist_agg_total, aes(x = year, y = disturbance_ha))

dist_agg_total_country <- dist_agg %>%
  group_by(year, country) %>%
  summarize(disturbance_ha = sum(disturbance_ha, na.rm = TRUE))

write_csv(dist_agg_total_country, "results/disturbance_area_map_based_total_countries.csv")

# Return intervals --------------------------------------------------------

mn <- mean(dist_agg_total[dist_agg_total$year < 2016, "disturbance_ha"][[1]])
sd <- sd(dist_agg_total[dist_agg_total$year < 2016, "disturbance_ha"][[1]])

prop <- 1 - pnorm(dist_agg_total$disturbance_ha, mn, sd)

prop_2018_2020 <- prod(prop[dist_agg_total$year >= 2018])

draws_binom <- rbinom(10000, 1, prop_2018_2020)

(length(draws_binom) + 1) / sum(draws_binom)

# Calculate disturbance anomaly ----------------------------------------------

dat <- dist_agg %>%
  filter(!is.na(gridindex)) %>%
  group_by(gridindex, year) %>%
  summarize(disturbance_ha = sum(disturbance_ha, na.rm = TRUE)) %>%
  ungroup() %>%
  split(.$gridindex) %>%
  map(~ right_join(., data.frame(gridindex = unique(.$gridindex),
                                 year = 1986:2020), 
                   by = c("gridindex", "year"))) %>%
  map(~ mutate_at(., .vars = vars(disturbance_ha), .funs = function(x) ifelse(is.na(x), 0, x))) %>%
  bind_rows()

dat <- dat %>%
  left_join(grid_eur %>% 
              st_drop_geometry() %>%
              group_by(gridindex) %>%
              summarise(country = COUNTRY[which.max(AREA_km2)],
                        iso_cc = ISO_CC[which.max(AREA_km2)]), by = c("gridindex"))

# Drop anomaly for Norway in 2020, as there is no data available
dat <- dat %>%
  filter(!(country == "Norway" & year == 2020))

forest <- forest_agg %>%
  group_by(gridindex) %>%
  summarise(forest_ha = sum(forest_ha),
            land_ha = sum(land_ha)) %>%
  mutate(forestcover = forest_ha / land_ha) %>%
  ungroup()

forest_grid <- grid %>%
  right_join(forest)

reference_period <- 1986:2015

dat <- dat %>%
  group_by(gridindex) %>%
  filter(sum(disturbance_ha) > 35) %>% # Exclude areas with less than 1 ha/yr of disturbances on average
  filter(sum(disturbance_ha[year %in% reference_period]) > 30) %>% # Exclude areas with less than 1 ha/yr of disturbances on average
  mutate(anomaly = disturbance_ha / mean(disturbance_ha[year %in% reference_period], na.rm = TRUE) - 1) %>% 
  ungroup()

dat <- dat %>%
  full_join(forest)

dat <- dat %>%
  mutate(forestcover = forest_ha / land_ha)

dat <- dat %>% filter(!is.na(forest_ha))

save(dat, file = "temp/dat.RData")

dat_grid <- grid %>%
  right_join(dat, by = "gridindex")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, st_crs(grid))
world <- st_crop(world, st_bbox(cntrs) + c(-0.05, -0.05, 0.01, 0.01) * as.double(st_bbox(cntrs)))

dat_grid_2018_2020 <- dat_grid %>%
  filter(year %in% 2018:2020) %>%
  mutate(disturbance_rate = disturbance_ha / forest_ha,
         disturbance_rate_capped = ifelse(disturbance_rate > 0.05, 0.05, disturbance_rate)) %>%
  # group_by(gridindex) %>%
  # summarise(anomaly = sum(anomaly),
  #           forestcover = unique(forestcover),
  #           forest_ha = unique(forest_ha)) %>%
  mutate(anomaly_capped = ifelse(anomaly > 5, 5, anomaly))

p_anomaly_map <- ggplot() +
  geom_sf(data = world, color = "black", fill = "lightgray") +
  geom_sf(data = dat_grid_2018_2020,
          aes(fill = anomaly_capped * 100), col = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  scale_fill_gradient2(low = "#2166ac", mid = "#FFFFFF", high = "#b2182b",
                       breaks = c(-100, 0, 100, 200, 300, 400, 500),
                       labels = c("-100%", "0%", "100%", "200%", "300%", "400%", ">500%")) +
  scale_color_gradient2(low = "#2166ac", mid = "#FFFFFF", high = "#b2182b",
                        breaks = c(-100, 0, 100, 200, 300, 400, 500),
                        labels = c("-100%", "0%", "100%", "200%", "300%", "400%", ">500%")) +
  theme_linedraw() +
  theme(panel.spacing = unit(0, "cm"),
        #panel.background = element_rect(fill = "#d1e5f0", color = "black", size = 1.125),
        panel.background = element_rect(fill = "white", color = "black", size = 1.125),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.125, "cm"),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 9, color = "black"),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(expand = FALSE, datum = NA) +
  labs(col = NULL, fill = NULL, 
       title = "a) Forest disturbance anomalies") +
  facet_wrap(~year, ncol = 3)

ggsave("results/disturbance_anomaly_cummulative_map_2018-2020.pdf", p_anomaly_map, width = 7.5, height = 3)

sum_temp_regions <- dat %>%
  filter(!is.na(disturbance_ha)) %>%
  left_join(cntr_groups, by = c("iso_cc" = "iso_code")) %>%
  filter(!is.na(euro_region)) %>%
  group_by(year, euro_region) %>%
  summarise(disturbance_ha = sum(disturbance_ha)) %>%
  group_by(euro_region) %>%
  mutate(anomaly = disturbance_ha / mean(disturbance_ha[year %in% reference_period], na.rm = TRUE) - 1) %>%
  ungroup() %>%
  mutate(euro_region = str_to_title(euro_region),
         euro_region_label = abbreviate(euro_region, 1))

p_anomaly_regions <- ggplot() +
  # geom_violin(data = sum_temp_regions,
  #             aes(x = reorder(euro_region_label, anomaly, mean),
  #                 y = anomaly * 100,
  #                 col = euro_region_label),
  #             alpha = 0.5, width = 0.2) +
  geom_jitter(data = sum_temp_regions,
              aes(x = reorder(euro_region_label, anomaly, mean), 
                  y = anomaly * 100, 
                  #col = euro_region_label,
                  fill = euro_region_label,
                  alpha = ifelse(year <= 2017, "a", "b"),
                  shape = ifelse(year <= 2017, "a", "b")),
              width = 0.15, stroke = 0, size = 2) +
  ggrepel::geom_text_repel(data = sum_temp_regions %>% filter(year > 2017),
                           aes(x = reorder(euro_region_label, anomaly, mean), 
                               y = anomaly * 100, 
                               label = year,
                               col = euro_region_label),
                           size = 2.5) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", size = 1.25),
        axis.text = element_text(size = 9, color = "grey30"),
        axis.title = element_text(size = 10)) +
  labs(x = NULL, y = "Forest disturbance anomaly (%)", col = NULL, fill = NULL) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  coord_flip() +
  scale_shape_manual(values = c(22, 21)) +
  scale_alpha_manual(values = c(0.2, 1))

ggsave("results/disturbance_anomaly_1986-2020_regions.pdf", p_anomaly_regions, width = 3.5, height = 3.5)

# Model -------------------------------------------------------------------

load(file = "temp/dat.RData")
climate_2018 <- read_csv("data/climate/era5_sm_vpd_summer_anomaly.csv")

modeldat <- dat %>%
  filter(!is.na(disturbance_ha)) %>% 
  filter(year %in% 2017:2020) %>%
  mutate(disturbance_rate = disturbance_ha / forest_ha,
         disturbance_count = as.integer(disturbance_ha / 0.09),
         forest_count = as.integer(forest_ha / 0.09)) %>%
  left_join(climate_2018 %>% filter(year %in% 2018:2020)) %>%
  left_join(climate_2018 %>% 
              filter(year %in% 2018) %>% 
              dplyr::select(gridindex, sm_z_18 = sm_z, vpd_z_18 = vpd_z)) %>%
  left_join(climate_2018 %>% 
              filter(year %in% 2019) %>% 
              dplyr::select(gridindex, sm_z_19 = sm_z, vpd_z_19 = vpd_z)) %>%
  left_join(climate_2018 %>% 
              filter(year %in% 2020) %>% 
              dplyr::select(gridindex, sm_z_20 = sm_z, vpd_z_20 = vpd_z)) %>%
  arrange(gridindex) %>%
  filter(!is.na(vpd_z)) %>% 
  filter(!is.na(sm_z))

save(modeldat, file = "temp/modeldat.RData")
load("temp/modeldat.RData")

modeldat_inp <- modeldat %>% 
  filter(disturbance_ha > 0) %>%
  mutate(year = factor(year))

fit0 <- lm(log(anomaly + 1) ~ sm_z_18 * vpd_z_18 * year, 
           data = modeldat_inp)

fit1 <- lm(log(anomaly + 1) ~ sm_z_18 * vpd_z * year, 
          data = modeldat_inp)

fit2 <- lm(log(anomaly + 1) ~ sm_z * vpd_z * year, 
          data = modeldat_inp)

AIC(fit0, fit1, fit2)

lmtest::lrtest(fit1, fit2)

summary(fit1)

broom::tidy(fit1) %>%
  write_csv(., "results/modelresults.csv")

prediction <- effects::effect("sm_z_18:vpd_z:year", 
                              fit1, 
                              xlevels = list(sm_z_18 = seq(-5, 5, length.out = 100),
                                             vpd_z = c(-1, 0, 1, 2))) %>%
  as.data.frame()

p_responsecurve <- ggplot(data = prediction) +
  geom_point(data = modeldat %>%
               filter(disturbance_ha > 0) %>%
               sample_frac(., 0.01),
               #mutate(vpd_z = cut(vpd_z, c(-1.5, -0.5, 0.5, 1.5, 2.5), labels = c(-1, 0, 1, 2))) %>%
               #filter(!is.na(vpd_z)),
             aes(x = sm_z_18, y = (anomaly) * 100),
             alpha = 0.1, shape = 1, fill = NA) +
  geom_ribbon(aes(x = sm_z_18, ymin = (exp(lower) - 1) * 100, ymax = (exp(upper) - 1) * 100, 
                  fill = factor(vpd_z)),
              alpha = 0.3) +
  geom_line(aes(x = sm_z_18, y = (exp(fit) - 1) * 100, col = factor(vpd_z))) +
  theme_classic() +
  scale_color_brewer(palette = "RdBu", direction = -1) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  theme(panel.background = element_rect(color = "black", size = 1.2),
        axis.text = element_text(size = 8, color = "grey30"),
        axis.title = element_text(size = 9),
        legend.position = "right",
        # legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9)) +
  labs(x = "Summer soil moisture anomaly in 2018", 
       y = "Disturbance anomaly (%)",
       col = "Summer vapor\npressure deficit\nanomaly", 
       fill = "Summer vapor\npressure deficit\nanomaly") +
  ylim(-100, 500) +
  facet_wrap(~year)
  #facet_grid(vpd_z~year)

ggsave("results/disturbance_anomaly_response_curve.pdf", p_responsecurve, width = 7.5, height = 2.5)

prediction <- effects::effect("sm_z_18:vpd_z_18:year", 
                              fit0, 
                              xlevels = list(sm_z_18 = seq(-4, 4, length.out = 100),
                                             vpd_z_18 = c(-1, 0, 1, 2))) %>%
  as.data.frame()

ggplot(data = prediction) +
  geom_point(data = modeldat %>%
               filter(disturbance_ha > 0) %>%
               sample_frac(., 0.01),
             aes(x = sm_z_18, y = (anomaly) * 100),
             alpha = 0.1, shape = 1, fill = NA) +
  geom_ribbon(aes(x = sm_z_18, ymin = (exp(lower) - 1) * 100, ymax = (exp(upper) - 1) * 100, 
                  fill = factor(vpd_z_18)),
              alpha = 0.1) +
  geom_line(aes(x = sm_z_18, y = (exp(fit) - 1) * 100, col = factor(vpd_z_18))) +
  theme_classic() +
  scale_color_brewer(palette = "RdBu", direction = -1) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  #facet_wrap(~year) +
  theme(panel.background = element_rect(color = "black", size = 1.2),
        axis.text = element_text(size = 8, color = "grey30"),
        axis.title = element_text(size = 9),
        legend.position = "right",
        # legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9)) +
  labs(x = "Summer soil moisture anomaly", 
       y = "Disturbance anomaly (%)",
       col = "Summer vapor\npressure anomaly", 
       fill = "Summer vapor\npressure anomaly") +
  ylim(-100, 500) +
  facet_wrap(~year)

prediction <- effects::effect("sm_z:vpd_z:year", 
                              fit2, 
                              xlevels = list(sm_z = seq(-4, 4, length.out = 100),
                                             vpd_z = c(-1, 0, 1, 2))) %>%
  as.data.frame()

ggplot(data = prediction) +
  geom_point(data = modeldat %>%
               filter(disturbance_ha > 0) %>%
               sample_frac(., 0.01),
             aes(x = sm_z, y = (anomaly) * 100),
             alpha = 0.1, shape = 1, fill = NA) +
  geom_ribbon(aes(x = sm_z, ymin = (exp(lower) - 1) * 100, ymax = (exp(upper) - 1) * 100, 
                  fill = factor(vpd_z)),
              alpha = 0.1) +
  geom_line(aes(x = sm_z, y = (exp(fit) - 1) * 100, col = factor(vpd_z))) +
  theme_classic() +
  scale_color_brewer(palette = "RdBu", direction = -1) +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  #facet_wrap(~year) +
  theme(panel.background = element_rect(color = "black", size = 1.2),
        axis.text = element_text(size = 8, color = "grey30"),
        axis.title = element_text(size = 9),
        legend.position = "right",
        # legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9)) +
  labs(x = "Summer soil moisture anomaly", 
       y = "Disturbance anomaly (%)",
       col = "Summer vapor\npressure anomaly", 
       fill = "Summer vapor\npressure anomaly") +
  ylim(-100, 500) +
  facet_wrap(~year)

### For supplement

dat %>%
  group_by(year, country) %>%
  summarise(disturbance_ha = sum(disturbance_ha)) %>%
  group_by(country) %>%
  mutate(disturbance_ha_reference1986to2915 = round(mean(disturbance_ha[year %in% reference_period]), 0),
         anomaly_percent = (disturbance_ha / disturbance_ha_reference1986to2915 - 1) * 100) %>%
  ungroup() %>%
  mutate(country = ifelse(country == "The Former Yugoslav Republic of Macedonia", "North Macedonia", country)) %>%
  filter(!is.na(country)) %>%
  filter(year %in% 2018:2020) %>%
  mutate(disturbance_ha_anomaly = paste0(disturbance_ha, ":", round(anomaly_percent, 2))) %>%
  dplyr::select(-disturbance_ha, -anomaly_percent) %>%
  spread(key = "year", value = "disturbance_ha_anomaly") %>%
  separate("2018", c("disturbance_ha_2018", "disturbance_anomaly_percent_2018"), "\\:") %>%
  separate("2019", c("disturbance_ha_2019", "disturbance_anomaly_percent_2019"), "\\:") %>%
  separate("2020", c("disturbance_ha_2020", "disturbance_anomaly_percent_2020"), "\\:") %>%
  write_excel_csv("results/disturbance_anomalies_map_estimates.csv")

# Climate data for Figure 1 -----------------------------------------------

climgrid <- read_sf("data/climate/climategrid_epsg3035.gpkg")
climgrid <- climgrid %>% 
  right_join(climate_2018)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, st_crs(cntrs))
world <- st_crop(world, st_bbox(cntrs) + c(-0.05, -0.05, 0.01, 0.01) * as.double(st_bbox(cntrs)))

sm_anomaly <- climgrid %>%
  filter(year %in% 2018:2020) %>%
  filter(gridindex %in% dat$gridindex) %>%
  dplyr::select(gridindex, vpd_z, sm_z, year) %>%
  mutate(sm_z = ifelse(sm_z < -4.5, -4.5, sm_z),
         sm_z = ifelse(sm_z > 4.5, -4.5, sm_z)) %>%
  ggplot(.) +
  geom_sf(data = world, color = "black", fill = "lightgray") +
  geom_sf(aes(fill = sm_z), col = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  scale_fill_gradient2(low = "#b2182b", mid = "#FFFFFF", high = "#2166ac",
                       breaks = c(-4.5, -3, -1.5, 0, 1.5, 3, 4.5),
                       limits = c(-4.5, 4.5),
                       labels = c("< -4.5", "-3", "-1.5", "0", "1.5", "3", "> 4.5")) +
  scale_color_gradient2(low = "#b2182b", mid = "#FFFFFF", high = "#2166ac",
                        breaks = c(-4.5, -3, -1.5, 0, 1.5, 3, 4.5),
                        limits = c(-4.5, 4.5),
                        labels = c("< -4.5", "-3", "-1.5", "0", "1.5", "3", "> 4.5")) +
  theme_linedraw() +
  theme(panel.spacing = unit(0, "cm"),
        #panel.background = element_rect(fill = "#d1e5f0", color = "black", size = 1.125),
        panel.background = element_rect(fill = "white", color = "black", size = 1.125),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.125, "cm"),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 9, color = "black"),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(expand = FALSE, datum = NA) +
  labs(col = NULL, fill = NULL,
       title = "b) Summer soil moisture anomalies") +
  facet_wrap(~year, ncol = 3)

ggsave("results/sm_anomaly.pdf", sm_anomaly, width = 7.5, height = 3)

 vpd_anomaly <- climgrid %>%
  filter(year %in% 2018:2020) %>%
  filter(gridindex %in% dat$gridindex) %>%
  dplyr::select(gridindex, vpd_z, sm_z, year) %>%
  mutate(vpd_z = ifelse(vpd_z > 4.5, 4.5, vpd_z)) %>%
  ggplot(.) +
  geom_sf(data = world, color = "black", fill = "lightgray") +
  geom_sf(aes(fill = vpd_z), col = NA) +
  geom_sf(data = world, color = "black", fill = NA) +
  scale_fill_gradient2(low = "#2166ac", mid = "#FFFFFF", high = "#b2182b",
                       breaks = c(-3, -1.5, 0, 1.5, 3, 4.5),
                       limits = c(-3, 4.5),
                       labels = c("-3", "-1.5", "0", "1.5", "3", "> 4.5")) +
  scale_color_gradient2(low = "#2166ac", mid = "#FFFFFF", high = "#b2182b",
                        breaks = c(-3, -1.5, 0, 1.5, 3, 4.5),
                        limits = c(-3, 4.5),
                        labels = c("-3", "-1.5", "0", "1.5", "3", "> 4.5")) +
  theme_linedraw() +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "white", color = "black", size = 1.125),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(0.125, "cm"),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 9, color = "black"),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(expand = FALSE, datum = NA) +
  labs(col = NULL, fill = NULL, 
       title = "c) Summer vapor pressure deficite anomalies") +
  facet_wrap(~year, ncol = 3)

ggsave("results/vpd_anomaly.pdf", vpd_anomaly, width = 7.5, height = 3)

p <- p_anomaly_map + sm_anomaly + vpd_anomaly + plot_layout(ncol = 1)

ggsave("results/figure01.pdf", p, width = 7.5, height = 8)
