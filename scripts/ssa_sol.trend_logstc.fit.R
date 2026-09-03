library(Rssa)
library(tidyverse)
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
trd_solar<-Ocean_solar_anomaly$trend.solar
L <- 64
s <- ssa(trd_solar, L = L)
# 1. Compute wcor
w <- Rssa::wcor(s, groups = 1:20)
# 2. Extract matrix and convert to long-format tibble
w_mat <- as.matrix(w)

w_df <- expand.grid(F1 = 1:20, F2 = 1:20) %>%
  mutate(value = abs(as.vector(w_mat)))

# 3. Plot with clear axes
ggplot(w_df, aes(x = F1, y = F2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "W-Correlation Matrix",
    x = "Eigenvector Index",
    y = "Eigenvector Index",
    fill = "W-Cor"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    panel.grid = element_blank()
  )
plot(s,type="values")
plot(s,type="paired",idx=1:10)


# 1. Extract the reconstructed trend component vector
rc <- reconstruct(s, groups = list(Trend = 1:2))
trend_vec <- as.numeric(rc$Trend)
Time = Ocean_solar_anomaly$dt.mnth
# 2. Build a data frame (replace time_vec with your actual time or index vector)
df_trend <- tibble(
  Time = Ocean_solar_anomaly$dt.mnth, # or your date/year axis: Ocean_solar_anomaly$time
  Trend = trend_vec,
  Original = trd_solar
)

# 3. Plot isolated trend overlaid on original input
ggplot(df_trend, aes(x = Time)) +
  geom_line(aes(y = Original), color = "black", linetype = 2, lwd = 1.2) +
  geom_line(aes(y = Trend), color = "firebrick", lwd = 0.6) +
  labs(
    title = "SSA Reconstructed Trend (Components 1–2)",
    x = "",
    y = "Anomaly / Trend "
  ) +
  theme_minimal()
#===========
# fit to generalized logistic
library(tidyverse)
library(minpack.lm)

# 1. Fit standard 4-parameter logistic model to the 1-2 SSA trend
fit_log <- nlsLM(
  Trend ~ L / (1 + exp(-k * (Time - t0))) + b,
  data = df_trend,
  start = list(L = 1.0, k = 0.03, t0 = 1990, b = -0.2)
)
round(coef(fit_log),3) #  L       k    t0        b
                    # 1.082  0.047 2015.870  -0.126
t0=coef(fit_log)[["t0"]]# 2015.870
# 2. Extract fitted industrial trend and compute residual
df_trend <- df_trend %>%
  mutate(
    industrial_logistic = predict(fit_log),
    natural_residual    = Trend - industrial_logistic
  )
which.max(diff(df_trend$industrial_logistic)) #1987
diff(df_trend$industrial_logistic)[1987]
# 3. Visualize the separation
ggplot(df_trend, aes(x = Time)) +
  geom_line(aes(y = Trend, color = "Trend.ssa"), size = 1) +
  geom_line(aes(y = industrial_logistic, color = "Industr."), size = 1, linetype = "dashed") +
  geom_line(aes(y = natural_residual, color = "Resid"), size = 0.9) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = t0,linetype = 2)+
  scale_color_manual(values = c(
    "Trend.ssa" = "firebrick",
    "Industr." = "black",
    "Resid" = "navy"
  )) +
  labs(
    title = "Separation :Industrial and Natural ",
    x = "",
    y = "Anomaly Units",
    color = "Comp.:"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("figs/separ_industr.natural.png")
