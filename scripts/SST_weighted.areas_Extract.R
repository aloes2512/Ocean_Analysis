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
# CROP BOTH the raw data:sst_stack.new and  active_weights.new
# --- Equatorial Atlantic Separation ---
atl_west <- crop(sst_stack.new, ext(309, 359, -5, 5))
atl_east <- crop(sst_stack.new, ext(0, 15, -5, 5))
atlantic_raster <- merge(atl_west, atl_east)
#-----Meditaranean Separation
# --- Mediterranean Separation ---
med_west <- crop(sst_stack.new, ext(353, 359.9, 30, 46))
med_east <- crop(sst_stack.new, ext(0, 36, 30, 46))
mediterranean_raster <- merge(med_west, med_east)
# "Equatorial_Indian"   = c(40, 100, -5, 5)
Equatorial_Indian_raster=crop(sst_stack.new,c(40, 100, -5, 5))
# NH SH Oceans
NH_raster=crop(sst_stack.new,c(0,359.9,0,90))
SH_raster=crop(sst_stack.new,c(0,359.9,-90,0))
