# data used solar locked anomalies
Ocean_sol_anomaly=readRDS("data/Ocean_solar_anom.rds")
summary(Ocean_sol_anomaly)
colnames(Ocean_sol_anomaly)
Ocean_sol_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Clean_Baseline_Trend))
#s_co2 <- ssa(co2, L = 84)
baseline.trnd=Ocean_sol_anomaly$Clean_Baseline_Trend
length(baseline.trnd) #2116
s_baseline=ssa(baseline.trnd,L=500)
#plot(s_ocean, type = "values")
plot(s_baseline,type = "values")
#W <- wcor(s_ocean,groups = 1:15)
W_baseline <- wcor(s_baseline,groups = 1:20)
as.vector(W_baseline)%>% length()
dfbase_wcor <- expand.grid(F1 = 1:20, F2 = 1:20)
dfbase_wcor$W_baseline <- as.vector(W_baseline)
#--------
ggplot(dfbase_wcor, aes(x = F1, y = F2, fill = W_baseline)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "W.Baseline-Correlation Matrix",
    x = "Component",
    y = "Component",
    fill = "|W-cor|"
  ) +
  theme_minimal(base_size = 10)
w_tidy <- as.data.frame(as.table(as.matrix(W_baseline)))
colnames(w_tidy) <- c("Comp1", "Comp2", "W_baseline")
w_tidy <- w_tidy %>%
  mutate(
    # Extract just the numbers from the component names (e.g., "F1" -> 1)
    Comp1 = parse_number(as.character(Comp1)),
    Comp2 = parse_number(as.character(Comp2)),
    Abs_Corr = abs(W_baseline)
  )
# 3. Create the ggplot
w_plot <- ggplot(w_tidy, aes(x = Comp1, y = Comp2, fill = Abs_Corr)) +
  geom_tile(color = "gray95", linewidth = 0.1) +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(1, 20, by = 5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(1, 20, by = 5), expand = c(0, 0)) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "W-Correlation Matrix",
    subtitle = "First 20 Components",
    x = "Component Index",
    y = "Component Index",
    fill = "|W-Corr|"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    # This rotates the horizontal axis text by 45 degrees and aligns them perfectly:
    axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )
w_plot
# reconstructed <- reconstruct(s_ocean, groups = list(Harmonics = 1:20))

# Plot the original difference against your clean SSA reconstruction
#plot(differnc, type = "l", col = "red", main = "SSA Reconstruction of Harmonics")
#lines(reconstructed$Harmonics, col = "blue", lwd = 2)
Reconstructd.trds=tibble(dt.mnth=Ocean_sol_anomaly$dt.mnth,
                    reconstructed = unlist(reconstruct(s_baseline,
                                  groups=list(Harmonics= 1:20))),
                    harmonics=      unlist(reconstruct(s_baseline,
                                  groups=list(Harmonics= 2:20))),
                    trend=          unlist(reconstruct(s_baseline,
                                  groups=list(trd=1))),
                    baseline.trnd=Ocean_sol_anomaly$Clean_Baseline_Trend)
Reconstructd.trds%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=baseline.trnd,col="sol.lckd.trd"))+
  geom_line(aes(y=harmonics,col="harmonics"))+
  geom_line(aes(y=trend,col="trend"))+
  labs(x="",title = "Decomposed Sol.-locked Trend")
ggsave("figs/Decomposed_sol.lckd.trd.tiff")
tibble(dt.mnth=Ocean_sol_anomaly$dt.mnth,
       baseline.trnd=Ocean_sol_anomaly$Clean_Baseline_Trend,
       extract_trd=unlist(reconstruct(s_baseline,groups= list(1))))%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=extract_trd,col="ssa.trnd"))+
  geom_line(aes(y=baseline.trnd,col="baseline"))+
  labs(x="",title = "Ocean Temp.Anomaly",
       subtitle = "solar locked anomalies decomposed\nby Singular Spectrum Analysis (ssa)")
ggsave("figs/decomposed.trd.tiff")
