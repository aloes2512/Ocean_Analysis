library(terra)
library(lubridate)
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
t=time(sst_stack.new) # length 2116
dt.mnth=year(t)+(month(t)-1)/12
sst_global=global(sst_stack.new*active_weights,"sum",na.rm=TRUE)[,1]
total_weight=global(active_weights,"sum",na.rm=T)%>%as.numeric()
sst_global=sst_global/total_weight
global_anomaly<-   tibble(dt.mnth=year(t)+(month(t)-1)/12,
                          anomaly= sst_global-mean(sst_global,na.rm=T))
# CROP BOTH the raw data:sst_stack.new and  active_weights.new
# --- Meditaranean Separation ---
# 1. Crop the two pieces of the Mediterranean from the native 0:359 stack
med_west_sst <- crop(sst_stack.new, ext(352, 359.9, 30, 46))
med_east_sst <- crop(sst_stack.new, ext(0, 36, 30, 46))

# 2. Crop the exact matching pieces from your global weights grid
med_west_wts <- crop(active_weights, ext(352, 359.9, 30, 46))
med_east_wts <- crop(active_weights, ext(0, 36, 30, 46))

# 3. Calculate the spatial sums for the time series (2116 rows)
sum_west <- global(med_west_sst * med_west_wts, "sum", na.rm = TRUE)[, 1]
sum_east <- global(med_east_sst * med_east_wts, "sum", na.rm = TRUE)[, 1]
#---total weights sum
total_wts_west <- sum(global(med_west_wts, "sum", na.rm = TRUE)[, 1])
total_wts_east <- sum(global(med_east_wts, "sum", na.rm = TRUE)[, 1])
total_area_weight <- total_wts_west + total_wts_east
mediterranean_time_series =tibble(t=t,
                                  medi_ts=(sum_west + sum_east) / total_area_weight)
# atlantic atl ext(-52, 16, -5, 5)) # adjusted to grab the full basin edge
atl_west_sst <- crop(sst_stack.new, ext(-52,359.9,-6,6))
atl_east_sst <- crop(sst_stack.new, ext(0, 16, -6, 6))
# 2. Crop the exact matching pieces from your global weights grid
atl_west_wts <- crop(active_weights, ext(-52,359.9,-6,6))
atl_east_wts <- crop(active_weights, ext(0, 16, -6, 6))
# sum atl west
sum_atl_west <- global(atl_west_sst * atl_west_wts, "sum", na.rm = TRUE)[, 1]
sum_atl_east <- global(atl_east_sst * atl_east_wts, "sum", na.rm = TRUE)[, 1]
total_wts_atl_west <- sum(global(atl_west_wts, "sum", na.rm = TRUE)[, 1])
total_wts_atl_east <- sum(global(atl_east_wts, "sum", na.rm = TRUE)[, 1])
total_atl_weight <- total_wts_atl_west + total_wts_atl_east
atlantic_time_series=tibble(t=t,
                            atl_ts=(sum_atl_west+sum_atl_east)/total_atl_weight)
# "Equatorial_Indian"   = c(40, 100, -5, 5)
indian_sst=crop(sst_stack.new,c(40, 100, -5, 5))
indian_wts <- crop(active_weights, ext(40, 100, -5, 5))
sum_indian <- global(indian_sst * indian_wts, "sum", na.rm = TRUE)[, 1]
total_wts_indian <- sum(global(indian_wts, "sum", na.rm = TRUE)[, 1])
indian_time_series=tibble(t=t,
                          indian_ts=sum_indian/total_wts_indian)
# NH  Oceans
NH_sst=crop(sst_stack.new,c(-1,359,0,90))
NH_wts <- crop(active_weights, ext(-1,359,0,90))
sum_NH= global(NH_sst*NH_wts,fun="sum",na.rm=TRUE)
total_wts_NH <- sum(global(NH_wts, "sum", na.rm = TRUE)[, 1])
nh_time_series=tibble(t=t,
                      nh_ts=sum_NH/total_wts_NH)
#SH Oceans
SH_sst=crop(sst_stack.new,c(-1,359,-90,0))
SH_wts <- crop(active_weights, ext(-1,359,-90,0))
sum_SH= global(SH_sst*SH_wts,fun="sum",na.rm=TRUE)
total_wts_SH <- sum(global(SH_wts, "sum", na.rm = TRUE)[, 1])
sh_time_series=tibble(t=t,
                      sh_ts=sum_SH/total_wts_SH)
# arctic ocean
arct_sst=crop(sst_stack.new,ext(-1,359,66,90))
arct_wts=crop(active_weights,ext(-1,359,66,90))
sum_arct=global(arct_sst*arct_wts,fun="sum",na.rm=TRUE)
total_wts_arct<- sum(global(arct_wts, "sum", na.rm = TRUE)[, 1])
# total_wts_arct 250.8427
arc_time_series=tibble(t=t,
                      arct_ts=sum_arct/total_wts_arct)
#antarctic
antarc_sst=crop(sst_stack.new,ext(-1,359,-89,-66))
antarc_wts=crop(active_weights,ext(-1,359,-89,-66))
sum_antarc=global(antarc_sst*antarc_wts,fun="sum",na.rm=TRUE)
total_wts_antarc<- sum(global(antarc_wts, "sum", na.rm = TRUE)[, 1])
# total weight antarctic 123.3528
antarc_time_series=tibble(t=t,
                          antarc_ts=sum_antarc/total_wts_antarc)

# add to exising areal ts
Area.ts=readRDS("data/Areas.ts.rds")
Area.ts$nh.ts=nh_time_series
Area.ts$sh.ts=sh_time_series
Area.ts$indian.ts=indian_time_series
Area.ts$medi.ts=mediterranean_time_series
Area.ts$atl.ts=atlantic_time_series
summary(Area.ts)
saveRDS(Area.ts,"data/Area8.ts.rds")
Area.ts=readRDS("data/Area8.ts.rds")
Area.ts$arct.ts=arc_time_series
Area.ts$antarc.ts=antarc_time_series
saveRDS(Area.ts,"data/Area10.ts.rds")
lst10=readRDS("data/Area10.ts.rds")
lst10$global.ts=global_anomaly
names(lst10)
saveRDS(lst10,"data/Area10+.ts.rds")
