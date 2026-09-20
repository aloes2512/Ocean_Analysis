library(tidyverse)
library(itsmr)
Ocean_NOAA.data<-readRDS("data/NOAA.OCEAN.ANOMALIES.rds")
o_df=Ocean_NOAA.data$data%>%dplyr::select(dt.mnth,"Anomaly"=anoma.mean)
dim(o_df) #2120
o_df<-o_df[1:2117,]
#  'obs_anomaly' is your 2120-month series
M.ssn=c("season",12,"season",6)
o_df<-o_df%>%mutate(anomaly=Resid(Anomaly,M.ssn))
obs_anomaly<-o_df$anomaly
dt.mnth=o_df$dt.mnth
# 1. Fit Polynomial (Degree 2)
fit_poly2 <- lm(obs_anomaly ~ dt.mnth+I(dt.mnth^2))

# 2.1 Fit LOESS (large span)
fit_loess <- loess(obs_anomaly ~ dt.mnth, span = 0.85)
# 2.2 Other fits
colnames(df_results)
# 3. Consolidate into a single comparison frame
df_comp <- data.frame(
  dt.mnth        = dt.mnth,
  Observed    = obs_anomaly,
  Richards    = df_results$richardsB, # From your 2-pass NLS fit
  Solar_Locked= df_results$nh_phase,     # Your TSI zero-crossing baseline
  LOESS_0.85  = predict(fit_loess),
  Poly_Deg2   = predict(fit_poly2)
) %>%
  pivot_longer(
    cols = -c(dt.mnth, Observed),
    names_to = "Method",
    values_to = "Trend_Value"
  )

# 4. Comparative Overlay Plot
ggplot() +
  geom_line(data = df_comp %>% subset(Method == "Richards"),
            aes(x = dt.mnth, y = Observed), color = "grey80", alpha = 0.5) +
  geom_line(data = df_comp,
            aes(x = dt.mnth, y = Trend_Value, color = Method, linetype = Method), linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Comparison of Secular Trend Extraction Methods",
    subtitle = "2117-Month Ocean Surface Temperature Anomaly Series",
    x = "Year",
    y = "Temperature Anomaly (°C)"
  ) +
  theme_minimal()
