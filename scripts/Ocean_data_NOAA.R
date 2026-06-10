library(terra)
library(tidyverse)
url.NOAA.psl="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"
#browseURL("https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/")
#download.file(url.NOAA.psl,
#              destfile = "data/sst.mnmean.nc",
#             method = "libcurl",
#             mode = "wb")
# download did overwrite sst.mnmean.nc
sst_stack <- rast("~/Downloads/sst.mnmean-2.nc")
dim(sst_stack) #89 180 2068
library(lubridate)
date=time(sst_stack) # "1854-01-01" "2026-04-01"
dt.mnth=year(date)+(month(date)-1)/12
length(dt.mnth) # 2068
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
Ocean_mean.anom=tibble(dt.mnth=time(sst_stack),
                       anoma.mean=ocean_anom_noaa )
# =========
# MEDIAN
library(matrixStats)
# Use the ORIGINAL sst_stack, not the multiplied one!
# 1. Combine ocean fraction and relative grid area into a single weight map
total_weights <- ocean_mask_na * weights_rel
# 2. Extract the weights as a simple 1D vector
w_vals <- terra::values(total_weights)[, 1]
# 3. Identify exactly which grid cells to keep (valid ocean cells)
# This will have a length of exactly 1658 based on your error message
keep_idx <- which(!is.na(w_vals) & w_vals > 0)
w_vals_clean <- w_vals[keep_idx]
sst_vals <- values(sst_stack, na.rm = FALSE)
full_sst_matrix <- terra::values(sst_stack, mat = TRUE)
# 5. Filter the matrix down to JUST your ocean cells
# This matrix will now have exactly 2068 rows
sst_matrix_clean <- full_sst_matrix[keep_idx, ]
# 6. Calculate the weighted median across time (columns)
# apply(..., 2) loops through each month. 'month_column' will have a length of 2068,
# perfectly matching 'w_vals_clean'.
global_medians <- apply(sst_matrix_clean, 2, function(month_column) {
  clean_idx <- !is.na(month_column)

  # If a month is entirely missing data (highly unlikely), return NA
  if (!any(clean_idx)) return(NA)

  matrixStats::weightedMedian(
    x = month_column[clean_idx],
    w = w_vals_clean[clean_idx]
  )
})


# 3. Extract the weights (ocean fraction)
ocean_mask_weighted=ocean_mask_na*weights_rel
# We flatten this to a single vector
# 2. Extract the weights as a simple 1D vector
w_vals <- terra::values(total_weights)[, 1]


is.na(w_vals)%>%sum() # 4714
which(!is.na(w_vals)) %>%length()# length 11306
which(w_vals>0)%>%unlist()%>%length() # 11306
# 3. Identify exactly which grid cells to keep (valid ocean cells)
# This will have a length of exactly 1658 based on your error message
keep_idx <- which(!is.na(w_vals) & w_vals > 0)

w_vals_clean <- w_vals[keep_idx]
length(w_vals_clean) # 11306
library(matrixStats)
# 4. Extract the entire SST stack into a standard R matrix
# Rows = All Spatial Grid Cells, Columns = Time (Months)
# This step is highly optimized in C++ and runs very fast
full_sst_matrix <- terra::values(sst_stack, mat = TRUE)# mat returns values as a matrix
# 5. Filter the matrix down to JUST your ocean cells
# This matrix will now have exactly 2068 columns
sst_matrix_clean <- full_sst_matrix[keep_idx, ]
dim(sst_matrix_clean)


# 2. Extract the corresponding weights
w_vector <- w_vals[keep_idx]
length(w_vector) # 11306
# 6. Calculate the weighted median across time (columns)
# apply(..., 2) loops through each month. 'month_column' will have a length of 2068,
# perfectly matching 'w_vals_clean'.
global_medians <- apply(sst_matrix_clean, 2, function(month_column) {
  clean_idx <- !is.na(month_column)

  # If a month is entirely missing data (highly unlikely), return NA
  if (!any(clean_idx)) return(NA)

  matrixStats::weightedMedian(
    x = month_column[clean_idx],
    w = w_vals_clean[clean_idx]
  )
})
global_medians%>%length()
# CALCULATE GLOBAL MEANS
# 1. Ensure any remaining runtime NAs (like dynamic sea ice) don't break the math.
# If your data has no NAs in the ocean cells, you can skip to step 2.
sst_matrix_zeroed <- sst_matrix_clean
sst_matrix_zeroed[is.na(sst_matrix_zeroed)] <- 0

# 2. Vectorized Weighted Mean using Matrix Multiplication (%*%)
# This multiplies the weights by the values and sums them for every month at once.
global_means <- as.vector((w_vals_clean %*% sst_matrix_zeroed) / sum(w_vals_clean))

NOAA.Ocean.anomalies=tibble(dt.mnth=dt.mnth,
                            date=time(sst_stack),
                       sst.mean=global_means,
                      sst.median=global_medians,
                      anoma.mean=sst.mean-mean(sst.mean,na.rm=T),
                      anoma.median=sst.median-mean(sst.median,na.rm=T))

summary(NOAA.Ocean.anomalies)




NOAA.Ocean.anomalies%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.median ),col=4)+
  geom_line(aes(y=anoma.mean ),col= 2)+

  labs(x="",y="anomaly K",title = "Ocean Temperature Anomalies",
       subtitle="Global Ocean Mean(blue) and Median(red)",caption = "data: noaa.ersst.v5.nc ")
summary(NOAA.Ocean.anomalies)
Ocean_temp.data=list(url.source="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc",
                     update="2026-05-16",
                     data.grid=sst_stack,
                     data=NOAA.Ocean.anomalies)

saveRDS(Ocean_temp.data,"data/NOAA.ocean.anomalies.rds")
Ocean_temp.data=readRDS("data/NOAA.ocean.anomalies.rds")
NOAA.Ocean.anomalies=Ocean_temp.data$data

