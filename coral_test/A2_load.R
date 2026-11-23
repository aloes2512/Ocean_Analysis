# A2_load.R : build annual, sign-corrected (up = warm) records -> out/records.rds
base <- "./"   # run with the working directory set to the coral_test folder
d <- paste0(base, "data/")
rd <- function(f, ...) read.table(paste0(d, f), header = TRUE, sep = "\t", comment.char = "#",
                                   na.strings = c("NA", "-999", "-99.9"), quote = "", fill = TRUE, ...)
ann <- function(t, x, minn = 1, yr = floor(t)) {   # annual mean by calendar year, minn obs required
  ok <- is.finite(x) & is.finite(t); yr <- yr[ok]; x <- x[ok]
  m <- tapply(x, yr, mean); n <- tapply(x, yr, length)
  m <- m[n >= minn]; data.frame(year = as.integer(names(m)), x = as.numeric(m))
}
R <- list()
add <- function(name, df, site, var, unit, sgn, yrconv, indep, note = "") {
  df <- df[is.finite(df$x), ]; df <- df[order(df$year), ]
  df$x <- df$x * sgn
  R[[name]] <<- list(name = name, df = df, site = site, var = var, unit = unit, sign = sgn,
                     yrconv = yrconv, indep = indep, note = note)
}
# Hon Tre monthly -> calendar-year mean
x <- rd("goodkin2021-tn-sst-month-noaa.txt"); a <- ann(x[[1]], x[[3]], 10)
add("HonTre_SST", a, "HonTre", "Sr/Ca SST (annual mean of monthly)", "degC", 1, "calendar year Jan-Dec", TRUE)
x <- rd("goodkin2021-tn-sst-seas-noaa.txt")
add("HonTre_wet", data.frame(year = x[[1]], x = x[[3]]), "HonTre", "Sr/Ca SST Sep-Nov", "degC", 1, "season Sep-Nov of year", FALSE, "same coral as HonTre_SST")
# Fiji
x <- rd("linsley2006-AB-noaa.txt")
add("Fiji_AB_d18O", data.frame(year = floor(x[[1]]), x = x[[2]]), "Fiji", "coral d18O (flipped)", "permil VPDB", -1, "annual, age_CE.5 = calendar year (Linsley 2006 chronology)", TRUE, "d18O includes salinity (SPCZ) signal")
x <- rd("linsley2006-1f-noaa.txt")
add("Fiji_1F_d18O", data.frame(year = floor(x[[1]]), x = x[[2]]), "Fiji", "coral d18O (flipped)", "permil", -1, "annual calendar", FALSE, "second core, same site")
add("Fiji_1F_SrCa_SST", data.frame(year = floor(x[[1]]), x = x[[6]]), "Fiji", "Sr/Ca SST", "degC", 1, "annual calendar", FALSE, "same core as 1F d18O")
x <- rd("fiji-tonga-idpo2008-noaa.txt")
add("FijiTonga_IDPO", data.frame(year = x[[1]], x = x[[2]]), "Fiji", "IDPO index from Fiji+Tonga d18O (5-core)", "dimensionless", 1, "annual", FALSE, "composite of Fiji/Tonga cores; dependent on Fiji_AB; sign set later vs Fiji_AB")
# Amedee New Caledonia
x <- rd("amedee2012_UpdateStephans2004_Anomalies-noaa.txt")
a <- ann(x$Year, x$MeanPAA, 9)
add("Amedee_SrCa", a, "NewCaledonia", "Sr/Ca anomaly (flipped: low Sr/Ca = warm)", "mmol/mol", -1, "calendar year of monthly values", TRUE)
# Rarotonga
x <- rd("linsley2006-2r-srca_noaa.txt"); a <- ann(x[[1]], x[[2]], 8)
add("Rarotonga_SrCa", a, "Rarotonga", "Sr/Ca (flipped)", "mmol/mol", -1, "calendar year of interpolated monthly values", TRUE)
x <- rd("linsley2006-2r-d18o_noaa.txt"); a <- ann(x[[1]], x[[2]], 8)
add("Rarotonga_d18O", a, "Rarotonga", "coral d18O (flipped)", "permil", -1, "calendar year", FALSE, "same core as Rarotonga_SrCa; salinity caveat")
# Dry Tortugas
x <- rd("delong2014-masterSST-noaa.txt"); a <- ann(x[[1]], x[[2]], 8)
add("DryTortugas_SST", a, "DryTortugas", "Sr/Ca SST, annual mean of monthly", "degC", 1, "calendar year", TRUE, "Atlantic (not Pacific)")
# Ryukyu
x <- rd("kawakubo2017-sst-annual.txt")
add("Ryukyu_SST", data.frame(year = x[[1]], x = x[[2]]), "Ryukyu", "Sr/Ca SST annual", "degC", 1, "annual (see source)", TRUE)
# Wilson 2010 Nino3.4
x <- rd("enso-wilson2010-comp-noaa.txt")
add("Wilson_N34_coral", data.frame(year = x[[1]], x = x[[2]]), "composite", "Nino3.4 reconstruction, corals only", "std units", 1, "annual", FALSE, "reuses corals: dependent")
add("Wilson_N34_tel", data.frame(year = x[[1]], x = x[[3]]), "composite", "Nino3.4 recon, teleconnected proxies (tree rings etc.)", "std units", 1, "annual", FALSE, "mostly non-coral, but not independent of ENSO-driven rings")
# Freund
x <- rd("freund2019enso.txt"); a <- ann(x$age_CE, x$NCT, 4, yr = x$age_CE)
add("Freund_NCT", a, "composite", "Nino Cold Tongue index, mean of 4 seasons", "index", 1, "age_CE with seasons 3,6,9,12 (MAM..DJF)", FALSE, "multi-proxy incl. corals")
# Tierney
x <- rd("tierney2015sst-wpacific-noaa.txt")
add("Tierney_WPac", data.frame(year = x[[1]], x = x[[2]]), "composite", "W Pacific coral SST composite", "degC anom", 1, "Apr-Mar", FALSE, "reuses corals")
x <- rd("tierney2015sst-epacific-noaa.txt")
add("Tierney_EPac", data.frame(year = x[[1]], x = x[[2]]), "composite", "E Pacific coral SST composite", "degC anom", 1, "Apr-Mar", FALSE, "reuses corals")
# comparison series
x <- read.table(paste0(base, "data/HadCET_meantemp_monthly_totals.txt"), skip = 5, header = FALSE, fill = TRUE)
x <- x[x[[1]] >= 1659 & x[[1]] <= 2025, ]
add("CET", data.frame(year = x[[1]], x = x[[14]]), "CET", "Central England annual mean T", "degC", 1, "calendar year", FALSE, "comparison")
x <- read.table(paste0(base, "data/pages2k_Full_ensemble_median_and_95pct_range.txt"), header = FALSE, sep = "\t", comment.char = "#", skip = 0, fill = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
x <- x[suppressWarnings(!is.na(as.numeric(x[[1]]))), ]; yy <- as.numeric(x[[1]]); v <- as.numeric(x[[3]])
add("PAGES2k_GMST", data.frame(year = yy[yy >= 1000 & yy <= 2000], x = v[yy >= 1000 & yy <= 2000]), "global", "PAGES2k ensemble median GMST", "degC", 1, "Apr-Mar", FALSE, "comparison")
g <- read.csv(paste0(base, "data/global_series.csv"))
a <- tapply(g$ersst[g$year <= 2025], g$year[g$year <= 2025], mean)
add("ERSSTv6_global", data.frame(year = as.integer(names(a)), x = as.numeric(a)), "global", "ERSSTv6 global mean SST anomaly", "degC", 1, "calendar year", FALSE, "comparison")
# sunspots
s <- read.table(paste0(base, "data/silso_monthly.csv"), sep = ";", header = FALSE)
names(s) <- c("y", "m", "t", "sn", "sd", "n", "def")
saveRDS(list(R = R, sun = s), paste0(base, "out/records.rds"))
for (r in R) cat(sprintf("%-18s %4d-%4d n=%d\n", r$name, min(r$df$year), max(r$df$year), nrow(r$df)))
