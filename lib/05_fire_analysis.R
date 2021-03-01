
# Packages ----------------------------------------------------------------

library(tidyverse)

# Load data ---------------------------------------------------------------

fire <- read_csv("data/fire/effis-ba-2019.csv")

fire <- fire %>%
  gather(key = country_iso, value = burntarea, -year)

regions <- read_csv("data/admin/country_grouping.csv")

disturbances <- read_csv("results/disturbance_area_map_based_total_countries.csv")

fire <- fire %>%
  left_join(regions, by = c("country_iso" = "iso_code"))

# Derive stats ------------------------------------------------------------

p <- fire %>%
  filter(euro_region %in% c("central", "north", "west", "south-east", "south-west", "east")) %>%
  mutate(country_name = ifelse(country_name == "The Former Yugoslav Republic of Macedonia", "North Macedonia", country_name)) %>%
  ggplot(., aes(x = year, y = burntarea)) +
  geom_line(col = "grey", size = 1) +
  facet_wrap(~country_name, scales = "free", ncol = 4) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", size = 1.2),
        axis.text = element_text(size = 8, color = "grey30"),
        axis.title = element_text(size = 9),
        strip.background = element_blank(),
        strip.text = element_text(size = 9)) +
  labs(x = "Year", y = "Burnt area (ha)")

ggsave("results/jrc_burnt_area.pdf", p, width = 7.5, height = 10)

p <- fire %>%
  filter(euro_region %in% c("central", "north", "west", "south-east", "south-west", "east")) %>%
  group_by(year, euro_region) %>%
  summarize(burntarea = sum(burntarea)) %>%
  mutate(euro_region = str_to_title(euro_region)) %>%
  ggplot(., aes(x = year, y = burntarea)) +
  geom_line(col = "grey", size = 1) +
  facet_wrap(~euro_region, scales = "free", ncol = 3) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", size = 1.2),
        axis.text = element_text(size = 8, color = "grey30"),
        axis.title = element_text(size = 9),
        strip.background = element_blank(),
        strip.text = element_text(size = 9)) +
  labs(x = "Year", y = "Burnt area (ha)")

ggsave("results/jrc_burnt_area_region.pdf", p, width = 7.5, height = 5)

fire %>%
  left_join(disturbances, by = c("country_name_short" = "country", "year")) %>%
  filter(year == 2018) %>%
  filter(euro_region %in% c("central", "north", "west", "south-east", "south-west", "east")) %>%
  mutate(percentage_burnt = (burntarea / disturbance_ha) * 100) %>%
  filter(euro_region == "west") %>%
  View(.)

fire %>%
  left_join(disturbances, by = c("country_name_short" = "country", "year")) %>%
  filter(year == 2018) %>%
  filter(euro_region %in% c("central", "north", "west", "south-east", "south-west", "east")) %>%
  mutate(percentage_burnt = (burntarea / disturbance_ha) * 100) %>%
  group_by(euro_region) %>%
  summarize(mean = mean(percentage_burnt),
            median = median(percentage_burnt),
            min = min(percentage_burnt),
            max = max(percentage_burnt))

fire %>%
  left_join(disturbances, by = c("country_name_short" = "country", "year")) %>%
  filter(year == 2018) %>%
  filter(euro_region %in% c("central", "north")) %>%
  summarize(burntarea = sum(burntarea),
            disturbance_ha = sum(disturbance_ha)) %>%
  mutate(percentage_burnt = (burntarea / disturbance_ha) * 100)

  
  
