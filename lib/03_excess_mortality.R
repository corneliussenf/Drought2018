
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

dist_agg_total_country <- dist_agg %>%
  group_by(year, country) %>%
  summarize(disturbance_ha = sum(disturbance_ha, na.rm = TRUE))

write_csv(dist_agg_total_country, "results/disturbance_area_map_based_total_countries.csv")

# Return intervals --------------------------------------------------------

mn <- mean(dist_agg_total[dist_agg_total$year < 2016, "disturbance_ha"][[1]])
sd <- sd(dist_agg_total[dist_agg_total$year < 2016, "disturbance_ha"][[1]])

prop <- 1 - pnorm(dist_agg_total$disturbance_ha, mn, sd)

prop[dist_agg_total$year >= 2018]

draws_binom_18 <- rbinom(10000, 1, prop[dist_agg_total$year == 2018])
draws_binom_19 <- rbinom(10000, 1, prop[dist_agg_total$year == 2019])
draws_binom_20 <- rbinom(10000, 1, prop[dist_agg_total$year == 2020])

(length(draws_binom_18) + 1) / sum(draws_binom_18)
(length(draws_binom_18) + 1) / sum(draws_binom_19)
(length(draws_binom_18) + 1) / sum(draws_binom_20)

# Calculate excess mortality ----------------------------------------------

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
  left_join(grid_eur %>% st_drop_geometry(), by = c("gridindex"))

# Drop anomaly for Norway in 2020, as there is no data avialbale
dat <- dat %>%
  filter(!(COUNTRY == "Norway" & year == 2020))

forest <- forest_agg %>%
  group_by(gridindex) %>%
  summarise(forest_ha = sum(forest_ha),
            land_ha = sum(land_ha)) %>%
  mutate(forestcover = forest_ha / land_ha) %>%
  ungroup()

forest_grid <- grid %>%
  right_join(forest)

ggplot() +
  geom_sf(data = forest_grid, aes(fill = forestcover > 0.40, col = forestcover > 0.40))

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

dat_grid <- grid %>%
  right_join(dat, by = "gridindex")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, st_crs(grid))
world <- st_crop(world, st_bbox(cntrs) + c(-0.05, -0.05, 0.01, 0.01) * as.double(st_bbox(cntrs)))

dat_grid_2018_2020 <- dat_grid %>%
  filter(year %in% 2018:2019) %>%
  group_by(gridindex) %>%
  summarise(anomaly = sum(anomaly),
            forestcover = unique(forestcover),
            forest_ha = unique(forest_ha)) %>%
  mutate(anomaly_capped = ifelse(anomaly > 6, 6, anomaly),
         anomaly_capped = ifelse(anomaly_capped < -4.5, -4.5, anomaly_capped))

p_anomaly_map <- ggplot() +
  geom_sf(data = world, color = "black", fill = "lightgray") +
  geom_sf(data = dat_grid_2018_2020, 
          aes(fill = anomaly_capped * 100, col = anomaly_capped * 100)) +
  geom_sf(data = world, color = "black", fill = NA) +
  scale_fill_gradient2(low = "#2166ac", mid = "#FFFFFF", high = "#b2182b",
                       breaks = c(-450, -300, -150, 0, 150, 300, 450, 600),
                       labels = c("-450%", "-300%", "-150%", "0%", "150%", "300%", "450%", ">600%")) +
  scale_color_gradient2(low = "#2166ac", mid = "#FFFFFF", high = "#b2182b",
                        breaks = c(-450, -300, -150, 0, 150, 300, 450, 600),
                        labels = c("-450%", "-300%", "-150%", "0%", "150%", "300%", "450%", ">600%")) +
  theme_linedraw() +
  theme(panel.spacing = unit(0, "cm"),
        panel.background = element_rect(fill = "#d1e5f0", color = "black", size = 1.25),
        legend.key.height = unit(1.75, "cm"),
        legend.key.width = unit(0.125, "cm"),
        # legend.position = "bottom",
        # legend.title.align = 0.5,
        # legend.direction = "horizontal",
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  # guides(fill = guide_colorbar(title = "Cummulative canopy mortality anomaly 2018-2020", title.position = "top"),
  #        color = guide_colorbar(title = "Cummulative canopy mortality anomaly 2018-2020", title.position = "top")) +
  coord_sf(expand = FALSE) +
  labs(col = NULL, fill = NULL, title = "Cummulative canopy mortality anomaly 2018-2020")

ggsave("results/disturbance_anomaly_cummulative_map_2018-2020.pdf", p_anomaly_map, width = 3.5, height = 3.5)

drought_threshold <- -2.5
non_drought_threshold <- -2.5

sum_temp_drought <- dat %>%
  left_join(climate_2018, by = c("gridindex")) %>%
  filter(!is.na(sm_z)) %>%
  group_by(year = year.x) %>%
  summarise(anomaly_mean = mean(anomaly),
            anomaly_sd = sd(anomaly),
            anomaly_se = anomaly_sd / sqrt(n()),
            anomaly_drought_mean = mean(anomaly[sm_z < drought_threshold]),
            anomaly_drought_sd = sd(anomaly[sm_z < drought_threshold]),
            anomaly_drought_se = anomaly_drought_sd / sqrt(sum(sm_z < drought_threshold)),
            anomaly_nondrought_mean = mean(anomaly[sm_z > non_drought_threshold]),
            anomaly_nondrought_sd = sd(anomaly[sm_z > non_drought_threshold]),
            anomaly_nondrought_se = anomaly_nondrought_sd / sqrt(sum(sm_z > non_drought_threshold))) %>%
  ungroup()

sum_temp_regions <- dat %>%
  left_join(cntr_groups, by = c("ISO_CC" = "iso_code")) %>%
  filter(!is.na(euro_region)) %>%
  group_by(year, euro_region) %>%
  summarise(anomaly_mean = mean(anomaly),
            anomaly_sd = sd(anomaly),
            anomaly_se = anomaly_sd / sqrt(n())) %>%
  ungroup() %>%
  mutate(euro_region = str_to_title(euro_region),
         euro_region_label = abbreviate(euro_region, 1))

sum_temp_regions_2018_2020 <- dat %>%
  left_join(cntr_groups, by = c("ISO_CC" = "iso_code")) %>%
  filter(!is.na(euro_region)) %>%
  filter(year >= 2018) %>%
  group_by(euro_region) %>%
  summarise(anomaly_mean = mean(anomaly),
            anomaly_sd = sd(anomaly),
            anomaly_se = anomaly_sd / sqrt(n())) %>%
  ungroup() %>%
  mutate(euro_region = str_to_title(euro_region),
         euro_region_label = abbreviate(euro_region, 1))

p_anomaly_drought <- ggplot() +
  geom_ribbon(data = sum_temp_drought,
              aes(x = year, 
                  ymin = (anomaly_mean - 2 * anomaly_se) * 100, 
                  ymax = (anomaly_mean + 2 * anomaly_se) * 100, 
                  fill = "All forest"), 
              alpha = 0.2) +
  geom_line(data = sum_temp_drought,
            aes(x = year, 
                y = anomaly_mean * 100, 
                col = "All forest"),
            alpha = 0.8) +
  geom_ribbon(data = sum_temp_drought,
              aes(x = year, 
                  ymin = (anomaly_drought_mean - 2 * anomaly_drought_se) * 100, 
                  ymax = (anomaly_drought_mean + 2 * anomaly_drought_se) * 100,
                  fill = "2018 drought area"), 
              alpha = 0.2) +
  geom_line(data = sum_temp_drought,
            aes(x = year, 
                y = anomaly_drought_mean* 100,
                col = "2018 drought area")) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey30") +
  theme_classic() +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.text = element_text(size = 9),
        panel.background = element_rect(color = "black", size = 1.25),
        axis.text = element_text(size = 9, color = "grey30"),
        axis.title = element_text(size = 10)) +
  labs(x = "Year", y = "Canopy mortality anomaly (%)", col = NULL, fill = NULL) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")

ggsave("results/disturbance_anomaly_1986-2020_drought.pdf", p_anomaly_drought, width = 3.5, height = 3.5)

p_anomaly_regions <- ggplot() +
  geom_boxplot(data = sum_temp_regions %>% filter(year <= 2015),
              aes(x = euro_region_label, y = anomaly_mean * 100, col = euro_region_label), alpha = 0.5, outlier.colour = NA, width = 0.2) +
  geom_jitter(data = sum_temp_regions %>% filter(year < 2018),
              aes(x = euro_region_label, y = anomaly_mean * 100, fill = euro_region_label),
              alpha = 0.25, shape = 21, stroke = 0, width = 0.15) +
  # geom_point(data = sum_temp_regions %>% filter(year >= 2018),
  #            aes(x = euro_region, y = anomaly_mean * 100, col = euro_region, fill = euro_region),
  #            size = 2, shape = 21, col = "black") +
  # ggrepel::geom_text_repel(data = sum_temp_regions %>% filter(year >= 2018),
  #                          aes(x = euro_region, y = anomaly_mean * 100, label = year),
  #                          size = 2) +
  geom_errorbar(data = sum_temp_regions_2018_2020,
                aes(x = euro_region_label, ymin = (anomaly_mean - anomaly_se * 2) * 100, ymax = (anomaly_mean + anomaly_se * 2) * 100, col = euro_region_label),
                width = 0) +
  geom_point(data = sum_temp_regions_2018_2020,
             aes(x = euro_region_label, y = anomaly_mean * 100, col = euro_region_label, fill = euro_region_label),
             size = 2, shape = 21, col = "black") +
  ggrepel::geom_text_repel(data = sum_temp_regions %>% 
                             filter(year <= 2015 & euro_region == "Central") %>% 
                             group_by(euro_region_label) %>% 
                             summarise(anomaly_mean = min(anomaly_mean)),
                           aes(x = euro_region_label, y = anomaly_mean * 100, label = "1986-2015"),
                           size = 2) +
  ggrepel::geom_text_repel(data = sum_temp_regions_2018_2020 %>% filter(euro_region == "Central"),
                           aes(x = euro_region_label, y = anomaly_mean * 100, label = "2018-2020"),
                           size = 2) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", size = 1.25),
        axis.text = element_text(size = 9, color = "grey30"),
        axis.title = element_text(size = 10)) +
  labs(x = NULL, y = "Canopy mortality anomaly (%)", col = NULL, fill = NULL) +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  coord_flip()

ggsave("results/disturbance_anomaly_1986-2020_regions.pdf", p_anomaly_regions, width = 3.5, height = 3.5)

# Stich plot together

layout <- "
AAAA#CCCCCCCC
BBBB#CCCCCCCC
"

p <- p_anomaly_drought + p_anomaly_regions + p_anomaly_map + plot_layout(design = layout)

ggsave("results/figure02.pdf", p, width = 7.5, height = 4.5)

### For supplement

sum_temp_cntr <- dat %>%
  group_by(year, country = COUNTRY) %>%
  summarise(anomaly_mean = mean(anomaly),
            anomaly_sd = sd(anomaly),
            anomaly_se = anomaly_sd / sqrt(n())) %>%
  ungroup() %>%
  mutate(country = ifelse(country == "The Former Yugoslav Republic of Macedonia", "North Macedonia", country))

sum_temp_cntr_2018_2020 <- dat %>%
  left_join(cntr_groups, by = c("ISO_CC" = "iso_code")) %>%
  filter(!is.na(euro_region)) %>%
  filter(year >= 2018) %>%
  group_by(country = COUNTRY) %>%
  summarise(anomaly_mean = mean(anomaly),
            anomaly_sd = sd(anomaly),
            anomaly_se = anomaly_sd / sqrt(n())) %>%
  ungroup() %>%
  mutate(country = ifelse(country == "The Former Yugoslav Republic of Macedonia", "North Macedonia", country))

p <- ggplot() +
  geom_jitter(data = sum_temp_cntr %>% filter(year < 2018),
              aes(x = country, y = anomaly_mean * 100, fill = country),
              alpha = 0.25, shape = 21, stroke = 0, width = 0.15) +
  geom_errorbar(data = sum_temp_cntr_2018_2020,
                aes(x = country, ymin = (anomaly_mean - anomaly_se * 2) * 100, ymax = (anomaly_mean + anomaly_se * 2) * 100, col = country),
                width = 0) +
  geom_point(data = sum_temp_cntr_2018_2020,
             aes(x = country, y = anomaly_mean * 100, col = country, fill = country),
             size = 2, shape = 21, col = "black") +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", size = 1.25),
        axis.text = element_text(size = 9, color = "grey30"),
        axis.title = element_text(size = 10)) +
  labs(x = NULL, y = "Canopy mortality anomaly (%)", col = NULL, fill = NULL) +
  coord_flip()

ggsave("results/hrv_country.pdf", p, width = 5.5, height = 5.5)

p <- ggplot() +
  geom_ribbon(data = sum_temp_cntr %>% filter(!is.na(country)),
              aes(x = year, 
                  ymin = (anomaly_mean - 2 * anomaly_se) * 100, 
                  ymax = (anomaly_mean + 2 * anomaly_se) * 100), 
              alpha = 0.2) +
  geom_line(data = sum_temp_cntr %>% filter(!is.na(country)),
            aes(x = year, 
                y = anomaly_mean * 100),
            alpha = 0.8) +
  facet_wrap(~country, scales = "free", ncol = 4) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", size = 1.2),
        axis.text = element_text(size = 8, color = "grey30"),
        axis.title = element_text(size = 9),
        strip.background = element_blank(),
        strip.text = element_text(size = 9)) +
  labs(x = "Year", y = "Canopy mortality anomaly (%)")

ggsave("results/disturbance_anomaly_1986-2020_countries.pdf", p, width = 7.5, height = 10)

p <- ggplot() +
  geom_ribbon(data = sum_temp_regions,
              aes(x = year,
                  ymin = (anomaly_mean - 2 * anomaly_se) * 100,
                  ymax = (anomaly_mean + 2 * anomaly_se) * 100,
              fill = euro_region,
              group = euro_region),
              alpha = 0.2) +
  geom_line(data = sum_temp_regions,
            aes(x = year,
                y = anomaly_mean * 100,
                col = euro_region,
                group = euro_region),
            alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey30") +
  theme_classic() +
  theme(legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.text = element_text(size = 9),
        panel.background = element_rect(color = "black", size = 1.25),
        strip.background = element_blank(),
        axis.text = element_text(size = 9, color = "grey30"),
        axis.title = element_text(size = 10)) +
  labs(x = "Year", y = "Canopy mortality anomaly (%)", col = NULL, fill = NULL) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~euro_region)

ggsave("results/disturbance_anomaly_1986-2020_regions_v2.pdf", p, width = 7.5, height = 5)
