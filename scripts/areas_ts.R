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
# CROP BOTH the raw data:sst_stack.new and  active_weights.new
# --- Equatorial Atlantic Separation ---
atl_west <- crop(sst_stack.new, ext(309, 359, -5, 5))
atl_east <- crop(sst_stack.new, ext(0, 15, -5, 5))
atlantic_raster <- merge(atl_west, atl_east)
weights_atl_west=crop(active_weights,ext(309, 359, -5, 5))
weights_atl_east=crop(active_weights,ext(0, 15, -5, 5))
weights_atl=merge(weights_atl_east,weights_atl_west)
weighted_Atlantic=atlantic_raster*weights_atl
#----------
atlantic_numerator <- global(weighted_Atlantic, fun = "sum", na.rm = TRUE)
atlantic_denominator <- global(weights_atl, fun = "sum", na.rm = TRUE)
atlantic.ts=tibble(t=time(sst_stack.new),
                 atlantic_ts = atlantic_numerator[, 1] / atlantic_denominator[1, 1]
)
#-----Meditaranean Separation
# --- Mediterean Separation ---
med_west <- crop(sst_stack.new, ext(353, 359.9, 30, 46))
med_east <- crop(sst_stack.new, ext(0, 36, 30, 46))
mediterranean_raster <- merge(med_west, med_east)
weights_med_west=crop(active_weights,ext(353, 359.9, 30, 46))
weights_med_east=crop(active_weights,ext(0, 36, 30, 46))
weights_med=merge(weights_med_west,weights_med_east)
weighted_Meditarian=mediterranean_raster*weights_med
meditaranian_numerator <- global(mediterranean_raster, fun = "sum", na.rm = TRUE)
meditaranian_denominator <- global(weights_med, fun = "sum", na.rm = TRUE)
meditaranian.ts=tibble(t=time(sst_stack.new),
                   meditaranian_ts = meditaranian_numerator[, 1] / meditaranian_denominator[1, 1]
)
# "Equatorial_Indian"   = c(40, 100, -5, 5)
Equatorial_Indian=crop(sst_stack.new,c(40, 100, -5, 5))
weights_Equatorial_Indian<- crop(active_weights, c(40, 100, -5, 5))
weighted_Equatorial_Indian<- Equatorial_Indian * weights_Equatorial_Indian
equatorial_indian_numerator=global(weighted_Equatorial_Indian,fun= "sum",na.rm=T)
equatorial_indian_denominator=global(weighted_Equatorial_Indian,fun="sum",na.rm=TRUE)
equatorial_indian.ts=tibble(t=time(sst_stack.new),
                            equatorial_indian.ts=equatorial_indian_numerator[,1]/weighted_Equatorial_Indian[1,1])
# NH SH Oceans
sst_nh_cropped=crop(sst_stack.new,c(0,359.9,0,90))
weights_nh_cropped <- crop(active_weights, c(0,359.9,0,90))
weighted_nh_stack <- sst_nh_cropped * weights_nh_cropped
sst_nh_numerator=global(weighted_nh_stack,fun="sum",na.rm=TRUE)
sst_nh_denominator=global(weighted_nh_stack,fun="sum",ba.rm=TRUE)
nh.ts=tibble(t=time(sst_stack.new),
             nh.ts=sst_nh_numerator/sst_nh_denominator[1,1])
#-----
SH_raster=crop(sst_stack.new,c(0,359.9,-90,0))
sst_sh_cropped=crop(sst_stack.new,c(0,359.9,-90,0))
weights_sh_cropped <- crop(active_weights, c(0,359.9,-90,0))
weighted_sh_stack <- sst_sh_cropped * weights_sh_cropped
sst_sh_numerator=global(weighted_sh_stack,fun="sum",na.rm=TRUE)
sst_sh_denominator=global(weights_sh_cropped ,fun="sum",na.rm=TRUE)
sh.ts=tibble(t=time(sst_stack.new),
             sh.ts=sst_sh_numerator/sst_sh_denominator[1,1])
Area.ts=readRDS("data/Areas.ts.rds")
Area.ts$nh.ts=nh.ts
Area.ts$sh.ts=sh.ts
Area.ts$equatorial_indian.ts=equatorial_indian.ts
Area.ts$meditaranian.ts=meditaranian.ts
Area.ts$atlantic.ts=atlantic.ts
summary(Area.ts)
saveRDS(Area.ts,"data/Area8.ts.rds")
