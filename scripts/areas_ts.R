library(terra)
# new data
sst_stack.new <- rast("data/sst.mnmean-3.nc")
berk_mask <- terra::rast("data/Land_and_Ocean_1x1_SST_at_sea_ice.nc")$land_mask
crs(sst_stack.new) <- crs(berk_mask)
# land fraction: mask has values between 0 and 1
ocean_frac <- 1 - berk_mask  # Ocean fraction [0,1]
# cell size depends on latitude min: 1741.084013, max: 49238.887519
weights_abs.new <- cellSize(sst_stack.new[[1]], unit="km")
weights_rel.new <- weights_abs.new / global(weights_abs.new, "max", na.rm=TRUE)[1,1]
ocean_mask_fixed.new <- terra::resample(ocean_frac, sst_stack.new, method = "near")

ocean_mask_na.new <- ocean_mask_fixed.new
# all complete land cells == NA:
ocean_mask_na.new[ocean_mask_na.new == 0] <- NA
# 4714 cells 2x2°
# --- THE CORRECTED WEIGHT LOGIC ---
# Define the combined spatial weight (Relative Area x Ocean Fraction)
# weights values correct cell size, ocean_mask fraction of land
total_weights.new <- weights_rel.new * ocean_mask_na.new
active_weights.new <- mask(total_weights.new, sst_stack.new[[1]])

# --- REGIONAL CROP TO NORTH ATLANTIC (AMO) ---
# correct  ext(-80, 0, 0, 65) to lon 0 : 359
# Definition matching  native -1 to 359 grid
ext_amo <- ext(280, 359, 0, 65)

# CROP BOTH the raw data:sst_stack.new and  active_weights.new
sst_amo_cropped <- crop(sst_stack.new, ext_amo)
weights_amo_cropped <- crop(active_weights.new, ext_amo)

# --- CONSERVATIVE REGIONAL MEAN CALCULATION ---
# Multiply the cropped data by the cropped weights
weighted_amo_stack <- sst_amo_cropped * weights_amo_cropped
# Calculate WEIGHTED MEAN
# Sum up the weighted values over time, and divide by the total active weight sum
amo_numerator <- global(weighted_amo_stack, fun = "sum", na.rm = TRUE)
amo_denominator <- global(weights_amo_cropped, fun = "sum", na.rm = TRUE)

# Your final pristine, area-weighted North Atlantic time series vector

AMO.ts=tibble(t=time(sst_stack.new),
       amo_ts = amo_numerator[, 1] / amo_denominator[1, 1]
)
# 2. North Pacific (PDO) Extent - cutting off the tropics
ext_pdo <- ext(110, 240, 20, 65)
sst_pdo_cropped <- crop(sst_stack.new, ext_pdo)
weights_pdo_cropped <- crop(active_weights.new, ext_pdo)
# Multiply the cropped data by the cropped weights
weighted_pdo_stack <- sst_pdo_cropped * weights_pdo_cropped
# Sum up the weighted values over time, and divide by the total active weight sum
pdo_numerator <- global(weighted_pdo_stack, fun = "sum", na.rm = TRUE)
pdo_denominator <- global(weights_pdo_cropped, fun = "sum", na.rm = TRUE)
PDO.ts=tibble(t=time(sst_stack.new),
              pdo_ts = pdo_numerator[, 1] / pdo_denominator[1, 1]
)
# 3. Tropical Pacific (ENSO) Extent
ext_nino34 <- ext(190, 240, -5, 5)
# Crop both the raw data and your true weight matrix
sst_nino34_cropped <- crop(sst_stack.new, ext_nino34)
weights_nino34_cropped <- crop(active_weights.new, ext_nino34)
# Multiply the cropped data by the cropped weights
weighted_nino34_stack <- sst_nino34_cropped * weights_nino34_cropped
# Sum up the weighted values over time, and divide by the total active weight sum
nino34_numerator <- global(weighted_nino34_stack, fun = "sum", na.rm = TRUE)
nino34_denominator <- global(weights_nino34_cropped, fun = "sum", na.rm = TRUE)
NINO34.ts=tibble(t=time(sst_stack.new),
              nino34_ts = nino34_numerator[, 1] / nino34_denominator[1, 1]
)
Areas.ts=list(amo=AMO.ts,
              pdo=PDO.ts,
              nino34=NINO34.ts)
saveRDS(Areas.ts,"data/Areas.ts.rds")
# add more areas
library(terra)

initialize_ocean_weights <- function(sst_path, mask_path) {
  message("Initializing SST stack and calculating spatial weights...")

  # 1. Load the raw datasets
  sst_stack <- rast(sst_path)
  berk_mask <- rast(mask_path)$land_mask

  # 2. Harmonize Coordinate Reference Systems
  crs(sst_stack) <- crs(berk_mask)

  # 3. Calculate ocean fraction from land mask
  ocean_frac <- 1 - berk_mask

  # 4. Calculate absolute and relative cell sizes based on latitude
  weights_abs <- cellSize(sst_stack[[1]], unit = "km")
  max_size    <- global(weights_abs, "max", na.rm = TRUE)[1, 1]
  weights_rel <- weights_abs / max_size

  # 5. Resample the mask to match the target SST resolution
  ocean_mask_fixed <- terra::resample(ocean_frac, sst_stack, method = "near")

  # 6. Set pure land cells (0% ocean) to NA to prevent them from skewing averages
  ocean_mask_na <- ocean_mask_fixed
  ocean_mask_na[ocean_mask_na == 0] <- NA

  # 7. Compute final combined spatial weights (Relative Area x Ocean Fraction)
  total_weights  <- weights_rel * ocean_mask_na
  active_weights <- mask(total_weights, sst_stack[[1]])

  # Return everything as a clean, named list
  return(list(
    sst_stack      = sst_stack,
    active_weights = active_weights
  ))
}
sst_path="data/sst.mnmean-3.nc"
mask_path="data/Land_and_Ocean_1x1_SST_at_sea_ice.nc"
sst_stack.new=initialize_ocean_weights(sst_path=sst_path,
                                       mask_path = mask_path)$sst_stack
active_weights=initialize_ocean_weights(sst_path=sst_path,
                                        mask_path = mask_path)$active_weights
#___________________

