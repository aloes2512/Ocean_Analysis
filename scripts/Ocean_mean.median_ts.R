# 1. Load packages and data
library(tidyverse)
library(lubridate)
library(terra)
library(here)
library(matrixStats)

Ocean_data_update <- readRDS(here("data", "Ocean_data_update.rds"))

download.time <- Ocean_data_update$downld.time
raw_path <- Ocean_data_update$raw_path
source.url <- Ocean_data_update$origin_url

sst_stack <- terra::unwrap(Ocean_data_update$sst_stack)
dim(sst_stack)

saveRDS(Ocean_data_update, here("data", "Ocean_data_update.rds"))

# 2. Build ocean fraction and area weight rasters
berk_mask <- terra::rast(here("data", "Land_and_Ocean_1x1_SST_at_sea_ice.nc"))$land_mask
dim(berk_mask)
crs(sst_stack) <- crs(berk_mask)

ocean_frac <- 1 - berk_mask
ocean_mask_fixed <- terra::resample(ocean_frac, sst_stack, method = "near")

ocean_mask_na <- ocean_mask_fixed
ocean_mask_na[ocean_mask_na == 0] <- NA

weights_abs <- terra::cellSize(sst_stack[[1]], unit = "km")
weights_rel <- weights_abs / terra::global(weights_abs, "max", na.rm = TRUE)[1, 1]
active_weights <- terra::mask(weights_rel, sst_stack[[1]])
total_weights <- ocean_mask_na * active_weights
ocean_mask_weighted <- total_weights

weighted_sst_stack <- sst_stack * total_weights
terra::writeRaster(weighted_sst_stack, here("data", "weighted_sst_stack.tif"), overwrite = TRUE)
terra::writeRaster(active_weights, here("data", "active_weights.tif"), overwrite = TRUE)

# 3. Weighted sum and weighted mean via terra::global()
sst.wsum <- terra::global(sst_stack, "sum", weights = total_weights, na.rm = TRUE)[, 1]
active_weight_sums <- terra::global(!is.na(sst_stack), "sum", weights = total_weights, na.rm = TRUE)[, 1]
sst.mean <- sst.wsum / active_weight_sums

# 4. Weighted median across time (columns)
w_vals <- terra::values(total_weights)[, 1]
keep_idx <- which(!is.na(w_vals) & w_vals > 0)
w_vals_clean <- w_vals[keep_idx]

full_sst_matrix <- terra::values(sst_stack, mat = TRUE)
sst_matrix_clean <- full_sst_matrix[keep_idx, , drop = FALSE]

global_medians <- apply(sst_matrix_clean, 2, function(month_column) {
  clean_idx <- !is.na(month_column)
  if (!any(clean_idx)) return(NA_real_)
  matrixStats::weightedMedian(
    x = month_column[clean_idx],
    w = w_vals_clean[clean_idx]
  )
})

# 5. Assemble anomalies tibble
date <- terra::time(sst_stack)
dt.mnth <- year(date) + (month(date) - 1) / 12

NOAA.Ocean.anomalies <- tibble(
  dt.mnth = dt.mnth,
  date = date,
  sst.wsum = sst.wsum,
  sst.mean = sst.mean,
  sst.median = global_medians,
  anoma.mean = sst.mean - mean(sst.mean, na.rm = TRUE),
  anoma.median = sst.median - mean(sst.median, na.rm = TRUE)
)

# 6. Plots
dir.create(here("figs"), showWarnings = FALSE, recursive = TRUE)

plt.glob.mean <- NOAA.Ocean.anomalies %>%
  ggplot(aes(x = dt.mnth)) +
  geom_line(aes(y = anoma.mean), col = "grey") +
  geom_point(aes(y = anoma.mean), size = 0.2, col = 2) +
  labs(
    x = NULL,
    y = "anomaly mean [K]",
    title = "Global Ocean Anomaly",
    caption = "psl.noaa.gov/Datasets/noaa.ersst.v6"
  )
ggsave(here("figs", "Global_Ocean_Anomaly.png"), plot = plt.glob.mean)

plt.anomaly_mean.median <- NOAA.Ocean.anomalies %>%
  ggplot(aes(x = dt.mnth)) +
  geom_line(aes(y = anoma.median), col = 4) +
  geom_line(aes(y = anoma.mean), col = 2) +
  labs(
    x = "",
    y = "anomaly K",
    title = "Ocean Temperature Anomalies",
    subtitle = "Global Ocean Mean (blue); Median (red)",
    caption = "data: noaa.ersst.v6.nc"
  )
ggsave(here("figs", "Global_Mean_Median.png"), plot = plt.anomaly_mean.median)

# 7. Save outputs and verify
Ocean_NOAA.data <- list(
  url.source = source.url,
  update = Sys.time(),
  data.grid = terra::wrap(sst_stack),
  data = NOAA.Ocean.anomalies
)

saveRDS(Ocean_NOAA.data, here("data", "NOAA.OCEAN.ANOMALIES.rds"))

saveRDS(
  list(
    global_stack = sst_stack,
    dates = date,
    grid = terra::wrap(sst_stack),
    ocean_mask_weighted = ocean_mask_weighted,
    Ocean_NOAA.data = Ocean_NOAA.data
  ),
  here("data", "processed_global_inputs.rds")
)

Ocean_NOAA.data_check <- readRDS(here("data", "NOAA.OCEAN.ANOMALIES.rds"))
NOAA.Ocean.anomalies_check <- Ocean_NOAA.data_check$data
sst_stack_check <- terra::unwrap(Ocean_NOAA.data_check$data.grid)
summary(Ocean_NOAA.data_check)
