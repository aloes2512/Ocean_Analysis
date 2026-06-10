# --- 00_ERSST_PreProcess.R ---
library(terra)
library(tidyverse)

# 1. LOAD RAW DATA
# Using terra for fast raster handling
sst_raw <- rast("data/ersst.v5.total.nc")
mask_berkeley <- rast("data/Berkeley_Land_Mask.nc") # 1x1 degree

# 2. ALIGN GRIDS
# Resample the 1x1 mask to the 2x2 ERSST grid using 'mean'
# to capture the fraction of ocean in each 2x2 cell.
mask_resampled <- resample(mask_berkeley, sst_raw, method = "bilinear")

# 3. AREA WEIGHTING
# Calculate cell area (km2) to account for latitudinal convergence
cell_area <- cellSize(sst_raw, unit = "km")

# 4. DEFINE AREAS OF INTEREST (AOI)
# Create a list of extents for your basins
aois <- list(
  Pacific_Eq  = ext(-180, -70, -10, 10),
  Indian_Eq   = ext(40, 100, -10, 10),
  Atlantic_Eq = ext(-50, 10, -10, 10),
  Baltic      = ext(10, 30, 53, 66)
)

# 5. EXTRACTION LOOP
# Apply mask, weight by area, and calculate the mean anomaly per time step
process_basin <- function(basin_ext) {
  cropped_sst  <- crop(sst_raw, basin_ext)
  cropped_mask <- crop(mask_resampled, basin_ext)
  cropped_area <- crop(cell_area, basin_ext)

  # The physical calculation: Sum(Value * Mask * Area) / Sum(Mask * Area)
  # This handles the "sophisticated" mask where coastal cells are partially ocean
  weighted_mean <- global(cropped_sst * cropped_mask * cropped_area, "sum", na.rm=TRUE) /
    global(cropped_mask * cropped_area, "sum", na.rm=TRUE)
  return(as.numeric(weighted_mean[,1]))
}

clean_data <- map(aois, process_basin)

# 6. SAVE INTERFACE
saveRDS(clean_data, "processed_data/basin_anomalies_list.rds")