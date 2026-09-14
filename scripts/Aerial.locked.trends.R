library(tidyverse)
library(itsmr)
library(splines)
library(lubridate)
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
data_list=readRDS("data/Area10.ts.rds")

regions=names(data_list)
solar_data=readRDS("data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))%>%dplyr::select(-TSI)
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
zero_crossing_times <- solar_data$dt.mnth[find_nodes(solar_data$SI)]
regions
df.global=analyze_ocean_region(data_list,
  region_col=regions[11],
  solar_data=solar_data,
  zero_crossing_times=zero_crossing_times,
  region_label = toupper(regions[11]))
n= length(data_list)

BaselineTrends<-vector("list",length=11)
names(BaselineTrends)=regions


BaselineTrends <- regions %>%
  purrr::map(~ {
    # .x represents the current region string loop iteration
    df <- analyze_ocean_region(data_list,
                               region_col = .x,
                               solar_data = solar_data,
                               zero_crossing_times = zero_crossing_times,
                               region_label = toupper(.x))

    df %>% dplyr::select(dt.mnth, Baseline_Trend, Region)
  }) %>%
  setNames(regions)


# Verify the structure of your northern hemisphere time series
summary(BaselineTrends$medi.ts)
BaselineTrends.long=bind_rows(BaselineTrends)
regions
colnames(BaselineTrends.long)
BaselineTrends.long<-BaselineTrends.long%>%
  mutate(Region=as_factor(Region))
BaselineTrends.long%>%subset(Region %in% c("NINO34","ATL.TS","INDIAN.TS"))%>%
  ggplot(aes(x=dt.mnth,y=Baseline_Trend,col = Region))+
  geom_line()+
  labs(x="",title="Compare Solar-Locked Trends",
       subtitle = "equatorial zones ± 5°")
ggsave("figs/Equatorial Solar Locked Trends.png")
# Meditaranean
BaselineTrends.long%>%subset(Region %in% c("MEDI.TS","AMO","PDO"))%>%
  ggplot(aes(x=dt.mnth,y=Baseline_Trend,col = Region))+
  geom_line()+
  labs(x="",title="Compare Midlatitude Trends",
       subtitle = "AMO, PDO, MEDI.TS")
ggsave("figs/Midlatitude Trends.png")
# add arctic and antarctic
Area.10.ts<-readRDS("data/Area10.ts.rds")
BaselineTrends.10<-vector("list",length=10)
names(BaselineTrends.10)<-c(names(BaselineTrends),"arct.ts","antarc.ts")
regions.10=names(BaselineTrends.10)
BaselineTrends.10 <- regions.10 %>%
  purrr::map(~ {
    # .x represents the current region string loop iteration
    df <- analyze_ocean_region(data_list=Area.10.ts,
                               region_col = .x,
                               solar_data = solar_data,
                               zero_crossing_times = zero_crossing_times,
                               region_label = toupper(.x))

    df %>% dplyr::select(dt.mnth, Baseline_Trend, Region)
  }) %>%
  setNames(regions)
# Aerial locked trends plot

Regions10_long.trend <- bind_rows(BaselineTrends.10)
Regions10_long.trend%>%subset(Region %in% c("ARCT.TS","ANTARC.TS"))%>%
  ggplot(aes(x=dt.mnth,y=Baseline_Trend,col=Region))+
  geom_line()+
  labs(x="",title = "Arctic- Antarctic Ocean Trend",
       subtitle= "area from polarcircles 66° polwards")


