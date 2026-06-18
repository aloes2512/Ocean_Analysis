library(dplyr)
library(lubridate)
library(itsmr)
library(splines)
library(ggplot2)

analyze_ocean_region <- function(data_list,
                                 region_col,
                                 solar_data,
                                 zero_crossing_times,
                                 region_label = "Region") {

  message(paste("Processing Baseline Trend for:", region_label))

  # 1. Extract the specific region's data frame from your list
  df_raw <- data_list[[region_col]]

  # 2. Extract the second column as a pure vector for itsmr processing
  ts_vector <- df_raw[, 2]%>%unlist()

  # 3. Replicate your precise time and anomaly transformations
  M <- c("season", 12, "season", 6)

  df_processed <- df_raw %>%
    mutate(
      # Center the raw time series using the isolated vector
      anoma_centered = ts_vector - mean(ts_vector, na.rm = TRUE),

      # Remove 12 and 6 month seasonality using itsmr
      anomaly        = Resid(anoma_centered, M),

      # Calculate continuous decimal months from the timestamp column 't'
      Year           = year(t),
      mnth           = month(t),
      dt.mnth        = Year + (mnth - 1) / 12
    ) %>%
    dplyr::select(dt.mnth, anomaly)

  # 4. Join with your monthly solar variation data
  df_combined <- df_processed %>%
    left_join(solar_data, by = "dt.mnth")

  # 5. Fit the continuous B-spline locked at your pre-calculated zero-crossings
  trend_fit <- lm(anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                  data = df_combined,
                  na.action = na.exclude)

  # 6. Extract the final solar-locked baseline trend
  df_combined <- df_combined %>%
    mutate(
      Baseline_Trend = predict(trend_fit),
      Region         = region_label
    )

  return(df_combined)
}
Areas.ts=readRDS("data/Areas.ts.rds")
data_list=Areas.ts
region_col="pdo"
solar_data=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))%>%dplyr::select(-TSI)
#zero_crossing_times
# find zero crossings of SI
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
nodes <- find_nodes(solar_data$SI)
# 1. Identify the time coordinates of your zero-crossings
zero_crossing_times <- solar_data$dt.mnth[find_nodes(solar_data$SI)]
region_label = "PDO"
# apply to
PDO_locked.trend=analyze_ocean_region(data_list=Areas.ts,
                     region_col="pdo",
                     solar_data = solar_data,
                     zero_crossing_times=zero_crossing_times,
                     region_label="PDO")
NINO34_locked.trend=analyze_ocean_region(data_list=Areas.ts,
                                                          region_col="nino34",
                                                          solar_data = solar_data,
                                                          zero_crossing_times=zero_crossing_times,
                                                          region_label="NINO34")
AMO_locked.trend=analyze_ocean_region(data_list=Areas.ts,
                                      region_col="amo",
                                      solar_data = solar_data,
                                      zero_crossing_times=zero_crossing_times,
                                      region_label="AMO")
PDO_locked.trend%>%ggplot(aes(x=dt.mnth,y=Baseline_Trend ))+geom_line()
# combine results
library(dplyr)
library(purrr)

# 1. Collect your outputs into a single list
regions_list <- list(
  AMO    = AMO_locked.trend,
  PDO    = PDO_locked.trend,
  NINO34 = NINO34_locked.trend
)

# 2. Combine them vertically into a single long row dataset
# This preserves the 'Region' column you built into the function!
Regions_long.trend <- bind_rows(regions_list)
Regions_long.trend%>%
  ggplot(aes(x=dt.mnth,y=Baseline_Trend,col=Region))+
  geom_line()
