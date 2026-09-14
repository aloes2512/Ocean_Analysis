# load the data
library(lubridate)
Ocean_data_update <- readRDS("data/Ocean_data_update.rds")
# Unwrap the raster object
download.time<-ymd(Ocean_data_update$downld.time)
source.url<-Ocean_data_update$raw_path
# data source
Ocean_data_update$file_url
sst_stack <- unwrap(Ocean_data_update$sst_stack)
dim(sst_stack)# 2120


# =========
#=============

berk_mask = terra::rast("data/Land_and_Ocean_1x1_SST_at_sea_ice.nc")$land_mask
crs(sst_stack)=crs(berk_mask)
# ocean cells to limit land included < 10%
ocean_mask <- terra::mask(berk_mask, berk_mask < 0.1,maskvalues=F)
ocean_frac <- 1 - berk_mask  # Ocean fraction [0,1]
weights_abs<- cellSize(sst_stack[[1]], unit="km")
weights_rel<- weights_abs / global(weights_abs, "max", na.rm=TRUE)[1,1]
ocean_mask_fixed<- terra::resample(ocean_frac, sst_stack, method = "near")

ocean_mask_na <- ocean_mask_fixed
ocean_mask_na[ocean_mask_na==0]<-NA
weighted_sst_stack<- sst_stack * ocean_mask_na*weights_rel
active_weights <- mask(weights_rel, sst_stack[[1]])

sum_of_ocean_weights <- global(active_weights, fun="sum", na.rm=TRUE)


ocean_mask_fixed <- terra::resample(ocean_frac, sst_stack, method = "near")
ocean_mask_na <- ocean_mask_fixed
ocean_mask_na[ocean_mask_na==0]<-NA
# MEDIAN
library(matrixStats)
# Use the ORIGINAL sst_stack, not the multiplied one!
# 1. Combine ocean fraction and relative grid area into a single weight map
total_weights <- ocean_mask_na* active_weights
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
ocean_mask_weighted=ocean_mask_na*total_weights
# We flatten this to a single vector
# 2. Extract the weights as a simple 1D vector
w_vals <- terra::values(weights_rel)[, 1]


is.na(w_vals)%>%sum() # 0
which(!is.na(w_vals)) %>%length()# length 16020
which(w_vals>0)%>%unlist()%>%length() # 16020
# 3. Identify exactly which grid cells to keep (valid ocean cells)
# This will have a length of exactly 1658 based on your error message
keep_idx <- which(!is.na(w_vals) & w_vals > 0)

w_vals_clean <- w_vals[keep_idx]
length(w_vals_clean) # 16020
library(matrixStats)
# 4. Extract the entire SST stack into a standard R matrix
full_sst_matrix <- terra::values(sst_stack, mat = TRUE)# mat returns values as a matrix
# 5. Filter the matrix down to JUST your ocean cells
sst_matrix_clean <- full_sst_matrix[keep_idx, ]
dim(sst_matrix_clean) # 16020 2120
# 2. Extract the corresponding weights
w_vector <- w_vals[keep_idx]
length(w_vector) # 16020
# 6. Calculate the weighted median across time (columns)
global_medians <- apply(sst_matrix_clean, 2, function(month_column) {
  clean_idx <- !is.na(month_column)
  # If a month is entirely missing data (highly unlikely), return NA
  if (!any(clean_idx)) return(NA)
 matrixStats::weightedMedian(
    x = month_column[clean_idx],
    w = w_vals_clean[clean_idx]
  )
})
global_medians%>%length() #2120
# CALCULATE GLOBAL MEANS
sst_matrix_zeroed <- sst_matrix_clean
sst_matrix_zeroed[is.na(sst_matrix_zeroed)] <- 0

# 2. Vectorized Weighted Mean using Matrix Multiplication (%*%)
global_means <- as.vector((w_vals_clean %*% sst_matrix_zeroed) / sum(w_vals_clean))
date=time(sst_stack) # 1850-01-01 2026-08-01

dt.mnth=year(date)+(month(date)-1)/12

NOAA.Ocean.anomalies=tibble(dt.mnth=dt.mnth,
                            date=time(sst_stack),
                            sst.mean=global_means,
                            sst.median=global_medians,
                            anoma.mean=sst.mean-mean(sst.mean,na.rm=T),
                            anoma.median=sst.median-mean(sst.median,na.rm=T))

NOAA.Ocean.anomalies%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.mean),col="grey")+
  geom_point(aes(y=anoma.mean),size=0.2,col=2)+
  labs(x=NULL,y="anomaly mean [K]",
       title = "Global Ocean Anomaly",
       caption = "psl.noaa.gov/Datasets/noaa.ersst.v5")
ggsave("figs/Global_Ocean_Anomaly.png")
#-----

NOAA.Ocean.anomalies%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anoma.median ),col=4)+
  geom_line(aes(y=anoma.mean ),col= 2)+

  labs(x="",y="anomaly K",title = "Ocean Temperature Anomalies",
       subtitle="Global Ocean Mean(blue) and Median(red)",caption = "data: noaa.ersst.v5.nc ")
summary(NOAA.Ocean.anomalies)
Ocean_NOAA.data=list(url.source="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc",
                     update="2026-09-13",
                     data.grid=wrap(sst_stack),
                     data=NOAA.Ocean.anomalies)

saveRDS(Ocean_NOAA.data,"data/NOAA.OCEAN.ANOMALIES.rds")
# check saving
Ocean_NOAA.data=readRDS("data/NOAA.OCEAN.ANOMALIES.rds")
NOAA.Ocean.anomalies=Ocean_NOAA.data$data
sst_stack<-unwrap(Ocean_NOAA.data$data.grid)

