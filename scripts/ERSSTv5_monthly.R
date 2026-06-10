# simplified template for a 200‑year oscillator
library(ncdf4)

mod_template <- function(t, A, f, phi, poly_coefs) {
  long_period <- A * sin(2 * pi * f * t + phi)
  p <- poly_coefs
  trend   <- p[1] + p[2]*t + p[3]*t^2 + p[4]*t^3 + p[5]*t^4
  trend + long_period
}
# Background of restored data
# merging Land Air Tamperature and interpolated SST
"https://doi.org/10.5281/zenodo.3634713"
browseURL("https://essd.copernicus.org/articles/12/3469/2020/#section3")
library(ncdf4)
#`gistemp1200_ghcnv4_ersstv5.nc
#is a netCDF climate dataset file containing the NASA GISS Surface Temperature Analysis (GISTEMP v4) global temperature anomalies. It merges meteorological station data (GHCNv4) and sea surface temperatures (ERSSTv5) at a
#smoothing radius.`
# Download ERSSTv5 or ERSSTv4 in R
# ERSSTv5 does not include land areas
nc.url="~/downloads/gistemp1200_GHCNv4_ERSSTv5.nc"
ERSST.nc=nc_open(nc.url)
anom_3d <- ncvar_get(ERSST.nc, "tempanomaly")
ts_mean   <- apply(anom_3d, 3, mean, na.rm = TRUE)
ref_date<-as.Date("1800-01-01")
time_raw  <- ncvar_get(ERSST.nc, "time") # time in days since "1800-01-01"
diff(time_raw)%>%mean() # average step distance 30.43675 days
dates     <- ref_date + as.integer(time_raw)   # adapt ref_date
range(dates)
ts_df <- data.frame(date = dates, anom = ts_mean)
ts_df <- ts_df[!is.na(ts_df$anom), ]
ERSST.nc
nc_close(ERSST.nc)
saveRDS(ts_df,"~/projects/Global_Temp/data/ERSSTv5.rds")
ts_df%>%ggplot(aes(x=date,y=anom))+geom_line()+
  labs(x="",title = "ERSSTv5 Anomaly")

ERSSTv5=list(data=ts_df,
             meta.dat=ERSST.nc)
ERSSTv5%>%saveRDS("~/projects/Global_Temp/ERSST.rds")
ERSSTv5$meta.dat%>%summary()
ERSSTv5$meta.dat[[8]]
ERSSTv5$meta.dat[[12]]
ERSSTv5$meta.dat[[15]]
# eliminate seasonal and linear trend
colnames(ts_df) #date anom
NROW(ts_df) # 1748
library(itsmr)
M <- c("season", 12, "trend", 1)
ts_df=ts_df%>%mutate(sst = Resid(anom, M))
colnames(ts_df) # date anom sst
ts_df%>%ggplot(aes(x=date,y=sst))+
  geom_line(col="grey")+geom_smooth(method = "loess",span= 0.5)
# combine with Berkley land mask
library(terra)
library(ncdf4)

# Load GISTEMP
berkley.path<- "~/Desktop/Klima_Energiewende/Daten/Berkley/"
list.files(berkley.path)[4]

gistemp <- rast("~/downloads/gistemp1200_GHCNv4_ERSSTv5.nc")  # 90×180 grid
res(gistemp)

# Load Berkeley 1° mask (land fraction 0-1)
berk_mask = rast("~/R-Studio Directory/zenodo_berkeley/Land_and_Ocean_1x1_SST_at_sea_ice.nc")$land_mask
# Extract layer
berk_mask <- aggregate(berk_mask, 2, mean)  # Resample to 2° to match GISTEMP
res(berk_mask)
# Ocean mask (land fraction ≤0.1)
ocean_mask <- berk_mask <= 0.1
# Or  ocean_mask= 1- berk_mask
ext(ocean_mask)
res(ocean_mask) # 1 1
ext(gistemp)

ocean_mask_1 <- 1-berk_mask
ocean_mask_2 <- aggregate(ocean_mask_1,2,mean)
# Apply
res(gistemp) # 2 2
res(ocean_mask_2) # 2 2
gistemp_ocean <- mask(gistemp, ocean_mask_2)

library(terra)

# Spatial mean for each time layer (returns data.frame)
ts_ocean <- global(gistemp_ocean, fun = "mean", na.rm = TRUE)

# Extract temperature values and dates
ocean_anom <- ts_ocean$mean  # Vector of global ocean means
dates <- time(gistemp_ocean) # POSIX dates from NetCDF
length(dates) # 1748
# Convert to time series object
ts_data <- ts(ocean_anom, start = 1880, frequency = 12)  # Monthly since ~1880
length(ts_data) #1748
giss_ocean.tbl=tibble(dates=dates,ts_data=ts_data)
giss_ocean.tbl%>%ggplot(aes(x=dates,y=ts_data))+geom_line()+
  labs(x="months", y= "SST anomalie",title = "gistemp1200_GHCNv4_ERSSTv5")
