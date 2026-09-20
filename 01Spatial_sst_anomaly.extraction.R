# summary compiled
# 01_Spatial data extraction
source("~/projects/Ocean_Analysis/scripts/Ocean_mean.median_ts.R")
ls(environment())
dim(NOAA.Ocean.anomalies) # 2120  6
rm(list = setdiff(ls(environment()), "NOAA.Ocean.anomalies"))
# data saved as list:
"data/NOAA.OCEAN.ANOMALIES.rds"
colnames(NOAA.Ocean.anomalies)#
#"dt.mnth"      "date"         "sst.mean"
#"sst.median"   "anoma.mean"   "anoma.median"
