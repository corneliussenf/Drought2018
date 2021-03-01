
# Libraries ---------------------------------------------------------------

library(raster)
library(tidyverse)
library(landscapemetrics)
library(foreach)
library(doParallel)
library(data.table)
library(rstudioapi)
library(sf)

rasterOptions(tmptime = 2, tmpdir = "temp/")
write("TMP = temp", file = file.path('~/.Renviron'))
source("lib/raster_as_data_table.R")
Rcpp::sourceCpp('lib/ngb_rcpp.cpp')

distmap_path <- "/home/csenf/Projects/mapping/results/version1.1"

# Settings ----------------------------------------------------------------

#cntr <- "albania"
countries <- list.files(distmap_path)
years <- 1986:2020
radius <- 5000

### Check if country was already processed, otherwise switch to next

cntr <- 0
k <- 0
while (cntr == 0) {
  k <- k + 1
  if (!file.exists(paste0("temp/patches/", countries[k]))) cntr <- countries[k]
}

print(cntr)

# Identify patches --------------------------------------------------------

print("Identifying final patches")

tictoc::tic()

dist_map <- raster(paste0(distmap_path, "/", cntr, "/disturbance_year_", cntr, ".tif"))

years_reclass <- raster::unique(dist_map)

inp_size <- file.size(paste0(distmap_path, "/", cntr, "/disturbance_year_", cntr, ".tif")) * 1e-6
cores <- ifelse(inp_size < 20, 20, ifelse(inp_size < 40, 10, 1))

if (cores == 1) {
  
  patches_centroids <- foreach(i = 1:length(years)) %do% {
    
    patch_dir <- paste0("temp/patches/", cntr)
    dir.create(patch_dir, recursive = TRUE, showWarnings = FALSE)
    
    dist_map_tmp <- reclassify(dist_map, matrix(c(years_reclass, years_reclass == years[i]), ncol = 2))
    
    patches <- clump(dist_map_tmp, gaps = FALSE, 
                     filename = paste0(patch_dir, "/patches_", years[i], "_", cntr, ".tif"), 
                     datatype = "INT4U", overwrite = TRUE)
    
    index_values <- expand.grid(col = 1:ncol(patches), row = 1:nrow(patches))
    index_values <- cbind(index_values, patch = getValues(patches))
    
    patches_centroids <- index_values %>%
      as.data.frame(.) %>%
      na.omit(.) %>%
      group_by(patch) %>%
      summarize(row = round(mean(row), 0),
                col = round(mean(col), 0)) %>%
      mutate(year = years[i])
    
    patches_centroids <- cbind(patches_centroids, xyFromCell(dist_map, cellFromRowCol(dist_map, patches_centroids$row, patches_centroids$col)))
    
    rm(patches)
    
    patches_centroids
    
  }
  
} else {
  
  cl <- makeCluster(cores)
  
  registerDoParallel(cl)
  
  patches_centroids <- foreach(i = 1:length(years), .packages = c("raster", "landscapemetrics", "tidyverse")) %dopar% {
    
    patch_dir <- paste0("temp/patches/", cntr)
    dir.create(patch_dir, recursive = TRUE, showWarnings = FALSE)
    
    dist_map_tmp <- reclassify(dist_map, matrix(c(years_reclass, years_reclass == years[i]), ncol = 2))
    
    patches <- clump(dist_map_tmp, gaps = FALSE, 
                     filename = paste0(patch_dir, "/patches_", years[i], "_", cntr, ".tif"), 
                     datatype = "INT4U", overwrite = TRUE)
    
    index_values <- expand.grid(col = 1:ncol(patches), row = 1:nrow(patches))
    index_values <- cbind(index_values, patch = getValues(patches))
    
    patches_centroids <- index_values %>%
      as.data.frame(.) %>%
      na.omit(.) %>%
      group_by(patch) %>%
      summarize(row = round(mean(row), 0),
                col = round(mean(col), 0)) %>%
      mutate(year = years[i])
    
    patches_centroids <- cbind(patches_centroids, xyFromCell(dist_map, cellFromRowCol(dist_map, patches_centroids$row, patches_centroids$col)))
    
    rm(patches)
    
    patches_centroids
    
  }
  
  stopCluster(cl)
  
}

patches_centroids <- patches_centroids %>%
  bind_rows() %>%
  as_tibble()

write_csv(patches_centroids, paste0("data/disturbances/patch_centroids/", cntr, "_patch_centroids.csv"))

gc()

print("Finished identifying final patches")

tictoc::toc()

# Calculate ngb metrics ---------------------------------------------------

print("Calculating neighborhood metrics")

tictoc::tic()

### Create kernel

maxd <- round(radius / 30)

dmat <- matrix(0, nrow= 2*maxd+1, ncol=2*maxd+1)

for (i in 1:(2 * maxd + 1))
  for (j in 1:(2 * maxd + 1))
    dmat[j, i] <- sqrt((maxd - i) * (maxd - i) + (maxd - j) * (maxd - j) ) * 30

dtab <- data.frame()

for (i in 1:(2 * maxd + 1))
  for (j in 1:(2 * maxd + 1))
    if (dmat[j, i] < radius)
      dtab <- rbind(dtab, data.frame(ix = j - maxd, iy = i - maxd))

### Loop through years and extract metrics

dist_map_matrix <- as.matrix(dist_map)

ngb_out <- vector("list", length(years))

for (i in 1:length(years)) {
  
  print(years[i])
  
  patches_centroids_tmp <- patches_centroids %>% filter(year == years[i])
  
  patches <- raster(paste0("temp/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))
  
  ngb <- countPxCentroidMT(dist_map_matrix, 
                           patches_centroids_tmp$col, 
                           patches_centroids_tmp$row,
                           patches_centroids_tmp$year,
                           dtab$ix, 
                           dtab$iy,
                           getValues(patches))
  
  colnames(ngb) <- c("n_t0", "n_total", "n_tminus1", "n_tplus1", "n_patches")
  ngb <- cbind(patches_centroids_tmp, ngb)
  ngb_out[[i]] <- as.data.frame(ngb)
  
}

### Write into final dataframe and save

ngb_metrics <- ngb_out %>%
  bind_rows() %>%
  mutate(perc_t0 = n_t0 / n_total,
         perc_tminus1 = n_tminus1 / n_total,
         perc_tplus1 = n_tplus1 / n_total) %>%
  as_tibble()

write_csv(ngb_metrics, paste0("data/disturbances/ngb/nbg_", cntr, ".csv"))

rm(patches)
rm(ngb)
rm(ngb_out)
rm(dist_map_matrix)
rm(dtab)
rm(dmat)

gc()

print("Finished calculating neighborhood metrics")

tictoc::toc()

# Patch metrics -----------------------------------------------------------

print("Calculating landscape metrics")

tictoc::tic()

### Calculate landscape metrics and save to disc

inp_size <- file.size(paste0(distmap_path, "/", cntr, "/disturbance_year_", cntr, ".tif")) * 1e-6
cores <- ifelse(inp_size < 20, 20, ifelse(inp_size < 40, 10, 1))

if (cores == 1) {
  
  
  ls_metrics <- foreach(i = 1:length(years)) %do% {
    
    patches <- raster(paste0("temp/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))
    edges <- boundaries(patches)
    
    lsm <- data.frame(patch = getValues(patches),
                      edges = getValues(edges)) %>%
      filter(!is.na(patch)) %>%
      group_by(patch) %>%
      summarize(area = n() * 900,
                perimeter = sum(edges) * 30) %>%
      mutate(frac = 2 * log(0.25 * perimeter) / log(area)) %>%
      mutate(year = years[i])
    
    gc()
    
    lsm
    
  }
  
} else {
  
  cl <- makeCluster(cores)
  
  registerDoParallel(cl)
  
  ls_metrics <- foreach(i = 1:length(years), .packages = c("raster", "tidyverse")) %dopar% {
    
    patches <- raster(paste0("temp/patches/", cntr, "/patches_", years[i], "_", cntr, ".tif"))
    edges <- boundaries(patches)
    
    lsm <- data.frame(patch = getValues(patches),
                      edges = getValues(edges)) %>%
      filter(!is.na(patch)) %>%
      group_by(patch) %>%
      summarize(area = n() * 900,
                perimeter = sum(edges) * 30) %>%
      mutate(frac = 2 * log(0.25 * perimeter) / log(area)) %>%
      mutate(year = years[i])
    
    gc()
    
    lsm
    
  }
  
  stopCluster(cl)
  
}

ls_metrics <- ls_metrics %>%
  bind_rows()

write_csv(ls_metrics, paste0("data/disturbances/lsm/lsm_", cntr, ".csv"))

gc()

print("Finished calculating landscape metrics")

tictoc::toc()

# Remove everything and restart R session ---------------------------------

rm(list = ls())

gc(reset = TRUE)

removeTmpFiles(h=0)

restartSession(command = 'source("lib/0X_spatialpattern.R")')



# Put everything together and plot results --------------------------------

names <- list.files("data/disturbances/ngb/") %>%
  gsub("nbg_", "", .) %>%
  gsub(".csv", "", .)

ngb <- list.files("data/disturbances/ngb/", full.names = TRUE) %>%
  map(read_csv) %>%
  #map(sample_frac, size = 0.01) %>%
  set_names(names) %>%
  bind_rows(., .id = "country")

names <- list.files("data/disturbances/lsm/") %>%
  gsub("lsm_", "", .) %>%
  gsub(".csv", "", .)

lsm <- list.files("data/disturbances/lsm/", full.names = TRUE) %>%
  map(read_csv) %>%
  #map(sample_frac, size = 0.01) %>%
  set_names(names) %>%
  bind_rows(., .id = "country")

# grid_eur <- read_sf("temp/grid_eur.gpkg")
# patches <- st_as_sf(ngb, coords = c("x", "y"))
# st_crs(patches) <- st_crs(grid_eur)
# patches_grid <- st_join(grid_eur, patches, join = st_intersects)
# patches_grid <- patches_grid %>% dplyr::select(gridindex, country, patch, year)
# save(patches_grid, file = "temp/patches_grid.RData")
load(file = "temp/patches_grid.RData")


# ngb <- ngb %>%
#   left_join(patches_grid)
 
lsm <- lsm %>%
  left_join(patches_grid)
 
# climate_2018 <- read_csv("data/climate/era5_sm_vpd_summer_anomaly.csv")

ngb_summary <- ngb %>%
  group_by(year) %>%
  summarise(ngb = mean(perc_t0))

# ngb_summary_drought <- ngb %>%
#   left_join(climate_2018) %>%
#   filter(sm_z < -2.5) %>%
#   group_by(year) %>%
#   summarise(ngb = mean(perc_t0))

# ngb_summary <- list(ngb_summary, ngb_summary_drought) %>%
#   set_names(c("allforest", "drought")) %>%
#   bind_rows(.id = "which")

forest_agg <- list.files("data/disturbances/aggregated_to_grid",
                         pattern = glob2rx("*forest*.csv$"), full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows()

freq_summary <- lsm %>%
  group_by(gridindex, year) %>%
  summarize(freq = n()) %>%
  ungroup() %>%
  left_join(forest_agg, by = "gridindex") %>%
  filter(forest_ha > 0) %>%
  mutate(freq = freq / forest_ha) %>%
  group_by(year) %>%
  summarize(freq = mean(freq))

lsm_summary <- lsm %>%
  group_by(year) %>%
  #summarise(lsm = mean(area / 10000))
  summarise(lsm = quantile(area / 10000, 0.95))

# lsm_summary_drought <- lsm %>%
#   left_join(climate_2018) %>%
#   filter(sm_z < -2.5) %>%
#   group_by(year) %>%
#   summarise(lsm = quantile(area / 10000, 0.95))
# 
# lsm_summary <- list(lsm_summary, lsm_summary_drought) %>%
#   set_names(c("allforest", "drought")) %>%
#   bind_rows(.id = "which")

dat <- ngb_summary %>% 
  left_join(lsm_summary) %>% 
  left_join(freq_summary)

save(dat, file = "temp/dat.RData")

library(plot3D)

pdf("results/ngb_lsm_freq.pdf", width = 5.5, height = 5.5)

scatter3D(x = dat$freq, y = dat$lsm, z = dat$ngb, colvar = dat$year, 
          phi = 35, 
          xlab = "Frequency", ylab = "Size", zlab = "Aggregation",
          type = "b", 
          lwd = 2,
          col = gg.col(), 
          colkey = FALSE)

text3D(x = dat$freq, y = dat$lsm, z = dat$ngb, colvar = dat$year, labels = ifelse(dat$year %in% c(1986, 2000, 2017:2020), dat$year, NA), 
       phi = 35, 
       xlab = "Frequency", ylab = "Size", zlab = "Aggregation",
       add = TRUE,
       col = gg.col(), 
       colkey = FALSE,
       cex = 0.75,
       adj = 1.25)

dev.off()

# For supplement

p <- dat %>%
  gather(key = key, value = value, -year) %>%
  mutate(key = case_when(
    key == "freq" ~ "Frequency (events per ha forest)",
    key == "lsm" ~ "Size (ha)",
    key == "ngb" ~ "Aggregation (0-1)"
  )) %>%
  split(.$key) %>%
  map(., ~ ggplot(., aes(x = year, y = value)) +
        geom_line(col = "grey", size = 1) +
        theme_classic() +
        theme(legend.position = "none",
              panel.background = element_rect(color = "black", size = 1.2),
              axis.text = element_text(size = 8, color = "grey30"),
              axis.title = element_text(size = 9),
              strip.background = element_blank(),
              strip.text = element_text(size = 9)) +
        labs(x = "Year", y = unique(.$key)))

p <- wrap_plots(p, ncol = 1)  

ggsave("results/resilience_indicators.pdf", p, width = 5.5, height = 7.5)

