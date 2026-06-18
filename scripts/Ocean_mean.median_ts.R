
library(terra)
library(tidyverse)
#sst_stack <- rast("~/Downloads/sst.mnmean-2.nc")
sst_stack.new<-rast("data/sst.mnmean-3.nc")
#dim(sst_stack.new) #89 180 2116
library(lubridate)
# =========
#=============
sst_stack.new<-rast("data/sst.mnmean-3.nc")
berk_mask = terra::rast("data/Land_and_Ocean_1x1_SST_at_sea_ice.nc")$land_mask
crs(sst_stack.new)=crs(berk_mask)
# ocean cells to limit land included < 10%
ocean_mask <- terra::mask(berk_mask, berk_mask < 0.1,maskvalues=F)
ocean_frac <- 1 - berk_mask  # Ocean fraction [0,1]
weights_abs.new <- cellSize(sst_stack.new[[1]], unit="km")
weights_rel.new <- weights_abs.new / global(weights_abs.new, "max", na.rm=TRUE)[1,1]
ocean_mask_fixed.new <- terra::resample(ocean_frac, sst_stack.new, method = "near")

ocean_mask_na.new <- ocean_mask_fixed.new
ocean_mask_na.new[ocean_mask_na.new==0]<-NA
weighted_sst_stack.new <- sst_stack.new * ocean_mask_na.new*weights_rel.new
active_weights.new <- mask(weights_rel.new, sst_stack.new[[1]])

sum_of_ocean_weights.new <- global(active_weights.new, fun="sum", na.rm=TRUE)


ocean_mask_fixed.new <- terra::resample(ocean_frac, sst_stack.new, method = "near")
ocean_mask_na.new <- ocean_mask_fixed.new
ocean_mask_na.new[ocean_mask_na.new==0]<-NA
# MEDIAN
library(matrixStats)
# Use the ORIGINAL sst_stack, not the multiplied one!
# 1. Combine ocean fraction and relative grid area into a single weight map
total_weights.new <- ocean_mask_na.new* active_weights.new
# 2. Extract the weights as a simple 1D vector
w_vals <- terra::values(total_weights.new)[, 1]
# 3. Identify exactly which grid cells to keep (valid ocean cells)
# This will have a length of exactly 1658 based on your error message
keep_idx <- which(!is.na(w_vals) & w_vals > 0)
w_vals_clean <- w_vals[keep_idx]
sst_vals <- values(sst_stack.new, na.rm = FALSE)
full_sst_matrix <- terra::values(sst_stack.new, mat = TRUE)
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
ocean_mask_weighted=ocean_mask_na.new*total_weights.new
# We flatten this to a single vector
# 2. Extract the weights as a simple 1D vector
w_vals <- terra::values(weights_rel.new)[, 1]


is.na(w_vals)%>%sum() # 5228
which(!is.na(w_vals)) %>%length()# length 10792
which(w_vals>0)%>%unlist()%>%length() # 10792
# 3. Identify exactly which grid cells to keep (valid ocean cells)
# This will have a length of exactly 1658 based on your error message
keep_idx <- which(!is.na(w_vals) & w_vals > 0)

w_vals_clean <- w_vals[keep_idx]
length(w_vals_clean) # 10792
library(matrixStats)
# 4. Extract the entire SST stack into a standard R matrix
# Rows = All Spatial Grid Cells, Columns = Time (Months)
# This step is highly optimized in C++ and runs very fast
full_sst_matrix <- terra::values(sst_stack.new, mat = TRUE)# mat returns values as a matrix
# 5. Filter the matrix down to JUST your ocean cells
# This matrix will now have exactly 2068 columns
sst_matrix_clean <- full_sst_matrix[keep_idx, ]
dim(sst_matrix_clean)


# 2. Extract the corresponding weights
w_vector <- w_vals[keep_idx]
length(w_vector) # 10792
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
date.new=time(sst_stack.new) # 1850-01-01 2026-04-01

dt.mnth.new=year(date.new)+(month(date.new)-1)/12

NOAA.Ocean.anomalies=tibble(dt.mnth=dt.mnth.new,
                            date=time(sst_stack.new),
                            sst.mean=global_means,
                            sst.median=global_medians,
                            anoma.mean=sst.mean-mean(sst.mean,na.rm=T),
                            anoma.median=sst.median-mean(sst.median,na.rm=T))

NOAA.Ocean.anomalies%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.mean),col="grey")+
  geom_point(aes(y=anoma.mean),size=0.2,col=2)
#-----

NOAA.Ocean.anomalies%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.median ),col=4)+
  geom_line(aes(y=anoma.mean ),col= 2)+

  labs(x="",y="anomaly K",title = "Ocean Temperature Anomalies",
       subtitle="Global Ocean Mean(blue) and Median(red)",caption = "data: noaa.ersst.v5.nc ")
summary(NOAA.Ocean.anomalies)
Ocean_NOAA.data=list(url.source="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc",
                     update="2026-05-16",
                     data.grid=sst_stack.new,
                     data=NOAA.Ocean.anomalies)

saveRDS(Ocean_NOAA.data,"data/NOAA.ocean.anomalies.2.rds")
rm(Ocean_NOAA.data)
Ocean_NOAA.data=readRDS("data/NOAA.ocean.anomalies.2.rds")
NOAA.Ocean.anomalies=Ocean_NOAA.data$data

