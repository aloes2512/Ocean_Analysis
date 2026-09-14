library(tidyverse)
library(gsignal)
library(splines)
library(itsmr)
Ocean_solar_anom<-readRDS("data/Ocean_solar_anom.rds")%>%
  dplyr::select(dt.mnth,Anomaly,"SI"=Solar_Variation,trd7)
Mssn=c("season",12,"season",6)
Ocean_solar_anom<-Ocean_solar_anom%>%mutate(anomaly=Resid(Anomaly,Mssn))
Ocean_solar_anom%>% ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=SI))
N<-NROW(Ocean_solar_anom)

find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
nodes <- find_nodes(Ocean_solar_anom$SI)
# time coordinates of zero-crossings
zero_crossing_times <- Ocean_solar_anom$dt.mnth[find_nodes(Ocean_solar_anom$SI)]

start_time <- min(Ocean_solar_anom$dt.mnth, na.rm = TRUE)
end_time   <- max(Ocean_solar_anom$dt.mnth, na.rm = TRUE)

# Keep interior knots strictly inside the dataset range
interior_knots <- zero_crossing_times[zero_crossing_times > start_time &
                                        zero_crossing_times < end_time]

fit_lm <- lm(
  anomaly ~ bs(dt.mnth,
               knots = interior_knots,
               Boundary.knots = c(start_time, end_time),
               degree = 3),
  data = Ocean_solar_anom
)
Ocean_solar_anom$trd.solar <- predict(fit_lm)
Ocean_solar_anom$secular_trend <- predict(
  loess(trd.solar ~ dt.mnth, data = Ocean_solar_anom, span = 0.5)
)
# ALTERATIVE
# Option A: span = 0.8 (141-year window)
Ocean_solar_anom$secular_trendA <- predict(
  loess(anomaly ~ dt.mnth, data = Ocean_solar_anom, span = 0.8)
)

# Option B: span = 1.0 (176-year window - recommended for preserving full Gleissberg wave)
Ocean_solar_anom$secular_trendB <- predict(
loess(anomaly ~ dt.mnth, data = Ocean_solar_anom, span = 1.0)
)
plt.secular<-Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=secular_trendA,col="A"))+
  geom_line(aes(y=secular_trendB,col="B"))+
  geom_line(aes(y=secular_trend,col="secular_trend"))+
  labs(x="",title = "Secular Trend ",
       subtitle = "secular_trend locked to irradiation\n A and B LOESS filtered Anomaly ")
ggsave("figs/LOESS & Solarlckd.png")
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd.solar,col="trd.solar"))+
  geom_line(aes(y=secular_trend,col="secular_trd"))
Ocean_solar_anom$absorbed_signal <- Ocean_solar_anom$anomaly - Ocean_solar_anom$secular_trend
Ocean_solar_anom<-Ocean_solar_anom%>%
  mutate(resids_solar=anomaly-trd.solar,
         resids_secular=anomaly-secular_trend)
Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=resids_solar),col=2)+
  geom_line(aes(y=resids_secular),col="grey")

Ocean_solar_anom%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=resids_solar),col="grey")+
  geom_line(aes(y=trd.solar,col="trd.solar"),linewidth = 1.3)+
  geom_line(aes(y=secular_trend,col="secular_trd"),linewidth = 1.3)

library(Rssa)
library(tidyverse)
# resids solar == difference of solar locked trend and anomaly
target_series<-Ocean_solar_anom$resids_secular
target_seriesB<-Ocean_solar_anom$anomaly-Ocean_solar_anom$secular_trendB
N <- length(target_series)
L_window <- 90*12 #  1080 month=90 years window

# Run SSA decomposition
s_ssa <- ssa(target_series, L = L_window, kind = "toeplitz-ssa")
s_ssa.B<-ssa(target_seriesB, L = L_window, kind = "toeplitz-ssa")


# Extract top 20 elementary components
n_comp <- 20
recs_all <- reconstruct(s_ssa, groups = as.list(1:n_comp))
recs_allB <- reconstruct(s_ssa.B, groups = as.list(1:n_comp))

# Find dominant period (in months) for each component using spec.pgram
comp_summary <- tibble(Component = 1:n_comp) %>%
  mutate(
    Variance_Pct = s_ssa$sigma[1:n_comp]^2 / sum(s_ssa$sigma^2) * 100,
    Dominant_Period_Months = map_dbl(1:n_comp, function(i) {
      s <- as.numeric(recs_all[[i]])
      spec <- spec.pgram(s, plot = FALSE)
      # Convert dominant frequency to period (in months)
      dom_freq <- spec$freq[which.max(spec$spec)]
      return(1 / dom_freq)
    })
  )
#-------
comp_summaryB <- tibble(Component = 1:n_comp) %>%
  mutate(
    Variance_Pct = s_ssa.B$sigma[1:n_comp]^2 / sum(s_ssa.B$sigma^2) * 100,
    Dominant_Period_Months = map_dbl(1:n_comp, function(i) {
      sB <- as.numeric(recs_allB[[i]])
      specB <- spec.pgram(sB, plot = FALSE)
      # Convert dominant frequency to period (in months)
      dom_freqB <- specB$freq[which.max(specB$spec)]
      return(1 / dom_freqB)
    })
  )
# Print the summary table to your R console
print(comp_summary, n = 20)


#========
#00000000
# 1. Prepare the target series (e.g., your hr() reconstructed wave or filtered residuals)
target_series <- Ocean_solar_anom$resids_secular # or resids_solar
solar_groups <- list(
  A= c(1,2),
  B= c(3,4),
  C = c(5, 6),
  D = c(7, 8)
)
solar_groupsB=list(
  a= c(1,2),
  b=c(3,4),
  c = c(5,6),
  d= c(7,8),
  e=c(9,12),
  f= c(10,11),
  g = c(13,14),
  h= c(17,19)
  )
# 2. Reconstruct components using Rssa
rec_solar <- reconstruct(s_ssa, groups = solar_groups)
rec_solarB <- reconstruct(s_ssa.B, groups = solar_groupsB)

# 3. Build a clean, tidy data frame for ggplot
df_rec <- Ocean_solar_anom %>%
  select(dt.mnth) %>%
  mutate(
    A = as.numeric(rec_solar$A),
    B = as.numeric(rec_solar$B),
    C = as.numeric(rec_solar$C),
    D = as.numeric(rec_solar$D),


  ) %>%
  pivot_longer(
    cols = -dt.mnth,
    names_to = "Mode",
    values_to = "Amplitude"
  )
colnames(df_rec)
df_recB <- Ocean_solar_anom %>%
  select(dt.mnth) %>%
  mutate(
    a = as.numeric(rec_solarB$a),
    b = as.numeric(rec_solarB$b),
    c = as.numeric(rec_solarB$c),
    d = as.numeric(rec_solarB$d),
    e= as.numeric(rec_solarB$e),
    f= as.numeric(rec_solarB$f),
    g = as.numeric(rec_solarB$g),
    h = as.numeric(rec_solarB$h),

  ) %>%
  pivot_longer(
    cols = -dt.mnth,
    names_to = "Mode",
    values_to = "Amplitude"
  )
# 4. Large, stacked ggplot with facets
ggplot(df_rec, aes(x = dt.mnth, y = Amplitude, color = Mode)) +
  geom_line(linewidth = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "SSA Decomposed Solar Cycles (Absorbed Signal)",
    x = "Date",
    y = "Temperature anomaly Contribution (°C)"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 16)
  )
#---------------
ggplot(df_recB, aes(x = dt.mnth, y = Amplitude, color = Mode)) +
  geom_line(linewidth = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "SSA Decomposed",subtitle="Secular Residuals",
    x = "",
    y = "anomaly contrib.(K)"
  )
#plot only selected
# A and B period 20 and 15 years
df_A.B <- Ocean_solar_anom %>%
  select(dt.mnth) %>%
  mutate(
    A = as.numeric(rec_solar$A),
    B = as.numeric(rec_solar$B),
    C = as.numeric(rec_solar$C),
    D = as.numeric(rec_solar$D)
  ) %>%
  pivot_longer(
    cols = -dt.mnth,
    names_to = "Mode",
    values_to = "Amplitude"
  )%>%subset(Mode%in% c("A","B") )
df_A.B<-df_rec%>%subset(Mode%in% c("A","B") )
df_A.B%>%ggplot(aes(x=dt.mnth,y = Amplitude, color = Mode))+
  geom_line(linewidth = 0.4)
df_rec%>%subset(Mode%in% c("C","D") )%>%
  ggplot(aes(x=dt.mnth,y = Amplitude, color = Mode))+
  geom_line(linewidth = 0.4)


#-----------
# plot correlation matrix
# Extract W-correlation matrix for the first 20 components
w_mat <- wcor(s_ssa, groups = 1:10)

# Convert matrix to a long data frame for ggplot
df_wcor <- expand.grid(F1 = 1:10, F2 = 1:10) %>%
  mutate(
    Correlation = as.vector(w_mat[1:10, 1:10]),
    F1 = factor(F1),
    F2 = factor(F2, levels = 1:10)   )

# Plot large heatmap
ggplot(df_wcor, aes(x = F1, y = F2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), size = 3) + # Print numerical correlation values
  scale_fill_gradient(low = "white", high = "firebrick") +
  theme_minimal(base_size = 14) +
  labs(
    title = "SSA W-Correlation Matrix (Groups 1 to 20)",
    x = "Component Number",
    y = "Component Number"
  ) +
  coord_fixed()


