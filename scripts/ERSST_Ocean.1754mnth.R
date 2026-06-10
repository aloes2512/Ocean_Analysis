fil.path="~/Desktop/Klima_Energiewende/Daten/gistemp250_GHCNv4.zarr/gistemp1200_GHCNv4_ERSSTv5.nc"
ERSST.nc=nc_open(fil.path)

#nc_close(ERSST.nc)
ersst= terra::rast(fil.path)
glob.temp.anom=tibble(dates=time(ersst),temp.anom=as.vector(global(ersst,mean,na.rm=T)[,1]))
summary(glob.temp.anom)
dates=time(ersst) # length 1754 month "1880-01-15" "2026-02-15"
# land/ ocean mask res 1,1 degree
berk_mask = terra::rast("~/R-Studio Directory/zenodo_berkeley/Land_and_Ocean_1x1_SST_at_sea_ice.nc")$land_mask
# reduce to res 2,2 to match ersst
berk_mask <- raster::aggregate(berk_mask, fact= 2, mean)  # 2
# reduce ocean mask to cells containing less than 10% land
ocean_mask <- terra::mask(berk_mask, berk_mask < 0.1,maskvalues=F)
crs(ocean_mask)==crs(ersst)
# get weighted average of ocean temp anomalies
ocean_frac <- 1 - berk_mask  # Ocean fraction [0,1]
# 1. Force the mask to match the SST grid perfectly
# 'near' is used because it's a categorical (0/1) mask
ocean_mask_fixed <- terra::resample(ocean_mask, ersst, method = "near")
plot(ocean_mask_fixed)
# 1. Get the area weights (accounts for cos(lat))
w <- cellSize(ersst, unit="km")
# 2. Normalize weights relative to the equator (optional but cleaner)
# This makes the equatorial cells = 1 and polar cells < 1
w_relative <- w / global(w, "max", na.rm=TRUE)[1,1]

final_weights <- w_relative * ocean_mask_fixed
weighted_ersst <- ersst * ocean_mask_fixed
# 1. Force the mask to be binary or NA
# Assuming Berkeley mask: 1 is ocean, 0 is land.
# We turn 0 into NA so it is completely ignored.
ocean_mask_strict <- ifel(ocean_mask_fixed > 0, 1, NA)

# 2. Apply the mask to your ERSST data
# This physically removes the land values from the computation
ersst_ocean_only <- mask(ersst, ocean_mask_strict)

# 3. Combine with weights for the Mean
# global() with na.rm=TRUE will now only look at non-NA ocean cells
final_mean <- global(ersst_ocean_only, fun="mean", weights=w_relative, na.rm=TRUE)
# Global mean ocean-weighted anomaly (e.g., for each time layer)
global_ocean.mean <- unlist(global(ersst_ocean_only, fun = "mean", na.rm = TRUE))
nc_close(ERSST.nc)
#======
Ocean_anoma=tibble(dates=time(ersst),ocean.anom=global_ocean.mean)
Ocean_anoma%>%ggplot(aes(x=dates,y=ocean.anom))+geom_line()
M.ssn=c("season",12)
Ocean_anoma=Ocean_anoma%>%mutate(yr.anom=Resid(ocean.anom,M.ssn),
                                 yr.smth=smooth.fft(yr.anom,0.1),
                                 res4=Resid(yr.smth,4),
                                 trd4=trend(yr.smth,4))
Ocean_anoma%>%ggplot(aes(x=dates,y=res4))+geom_line()
