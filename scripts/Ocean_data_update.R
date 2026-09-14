# update Ocean Mean Median data
#url.noaa<-"https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"
#browseURL("https://psl.noaa.gov/data/gridded/index.html")

dest_file <- "data/sst.mnmean-v6.nc"
# Direct URL to the single aggregated monthly file
file_url <- "https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v6/sst.mnmean.nc"
# Download with increased timeout for larger files
options(timeout = 600)
download.file(
  url = file_url,
  destfile = dest_file,
  mode = "wb" # 'wb' binary mode is required for NetCDF on Windows/macOS
)
#------
library(terra)
library(tidyverse)

sst_stack <- rast("data/sst.mnmean-v6.nc")
dim(sst_stack) #89 180 2120 download 20260913
Ocean_data_update= list(raw_path="~/projects/Ocean_Analysis/data/sst.mnmean-v6.nc",
                        downld.time="2026-09-13",
                        sst_stack=wrap(sst_stack),
                        file_url ="https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v6/sst.mnmean.nc"
)
saveRDS(Ocean_data_update,"~/projects/Ocean_Analysis/data/Ocean_data_update.rds")
# check
# Load the RDS file
Ocean_data_update <- readRDS("data/Ocean_data_update.rds")

# Unwrap the raster object
sst_stack <- unwrap(Ocean_data_update$sst_stack)

# Verify it works
dim(sst_stack)
#-----

