library(terra)
library(tidyverse)
url.NOAA.psl="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"

download.file(url.NOAA.psl,
              destfile = "data/sst.mnmean.nc",
              method = "libcurl",
              mode = "wb")

sst_stack <- rast("~/Downloads/sst.mnmean-2.nc")
#dim(sst_stack) #89 180 2068
dt.mnth=time(sst_stack) # "1854-01-01" "2026-04-01"
#length(dt.mnth) # 2068
berk_mask = terra::rast("data/Land_and_Ocean_1x1_SST_at_sea_ice.nc")$land_mask
crs(berk_mask)=crs(sst_stack)
# ocean cells to limit land included < 10%
ocean_mask <- terra::mask(berk_mask, berk_mask < 0.1,maskvalues=F)
ocean_frac <- 1 - berk_mask  # Ocean fraction [0,1]
# 1. Force the mask to match the SST grid perfectly
# 'near' is used because it's a categorical (0/1) mask
ocean_mask_fixed <- terra::resample(ocean_frac, sst_stack, method = "near")
# Convert 0 to NA so they are excluded from the "ranking"
ocean_mask_na <- ocean_mask_fixed
ocean_mask_na[ocean_mask_na == 0] <- NA
#________
weights_abs <- cellSize(sst_stack[[1]], unit="km")
# 2. Convert to relative weights (0 to 1 scale)
# This divides every cell by the area of an equatorial cell
weights_rel <- weights_abs / global(weights_abs, "max", na.rm=TRUE)[1,1]
#_________
weighted_sst_stack <- sst_stack * ocean_mask_na*weights_rel
#=====
# not global mean corrected mean taking only ocean cells
# Correct Mean = Sum of weighted values / Sum of weights
#ocean.anom_noaa<- unlist(global(weighted_sst_stack, fun="mean", na.rm=TRUE))
ocean.anom_noaa.sum=unlist(global(weighted_sst_stack, fun="sum", na.rm=TRUE))
# Create a weight layer that only exists where the SST data exists
active_weights <- mask(weights_rel, sst_stack[[1]])
# Now your denominator is perfectly synced to the available data
sum_of_ocean_weights <- global(active_weights, fun="sum", na.rm=TRUE)
ocean_anom_noaa <- ocean.anom_noaa.sum / as.numeric(sum_of_ocean_weights)
# =========
# MEDIAN
library(matrixStats)
# Use the ORIGINAL sst_stack, not the multiplied one!
sst_vals <- values(sst_stack, na.rm = FALSE)
# 3. Extract the weights (ocean fraction)
ocean_mask_weighted=ocean_mask_na*weights_rel
# We flatten this to a single vector
w_vals <- values(ocean_mask_weighted, na.rm = FALSE)
is.na(w_vals)%>%sum() # 4714
which(!is.na(w_vals))
which(w_vals>0)%>%unlist()%>%length() # 11306
# 4. Clean the data (Crucial Step!)
# Remove cells that are NA in either the temperature or the mask
# Also remove cells where ocean fraction is 0 (land)
keep_idx <- which(!is.na(w_vals) & w_vals > 0)
# 5. Apply the weighted median across each layer (column)
library(matrixStats)

# 1. Extract the values for just your 4,714 ocean cells
# This creates a matrix where:
# Rows = 4,714 space points
# Columns = ~2,068 time points (months)
sst_matrix <- sst_stack[keep_idx]

# 2. Extract the corresponding weights
w_vector <- w_vals[keep_idx]

# 3. Calculate the Weighted Median for every column (month) at once
# matrixStats::colWeightedMedians is optimized for this exact task
ocean_medians <- colWeightedMedians(
  as.matrix(sst_matrix),
  w = w_vector,
  na.rm = TRUE
)

# 4. Attach the dates for your physical model
results_df <- data.frame(
  date = time(sst_stack),
  median_sst = ocean_medians
)

ocean.median.anom=ocean_medians- mean(ocean_medians)

NOAA.Ocean.anomalies=tibble(dt.mnth=time(sst_stack),
                            mean.anom=unlist(ocean_anom_noaa)-mean(unlist(ocean_anom_noaa)) ,
                            median.anom=ocean.median.anom)
NOAA.Ocean.anomalies%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=median.anom),col=4)+
  geom_line(aes(y=mean.anom),col= 2)+

  labs(x="",y="anomaly K",title = "Ocean Temperature Anomalies",
       subtitle="Global Ocean Mean(red) and Median(blue)",caption = "data: noaa.ersst.v5.nc ")
summary(NOAA.Ocean.anomalies)
Ocean_temp.data=list(url.source="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc",
                     update="2026-05-16",
                     data.grid=sst_stack,
                     data=NOAA.Ocean.anomalies)

saveRDS(Ocean_temp.data,"data/NOAA.ocean.anomalies.rds")
Ocean_temp.data=readRDS("data/NOAA.ocean.anomalies.rds")
NOAA.Ocean.anomalies=Ocean_temp.data$data
##
