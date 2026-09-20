# ssa trends
head(oc)
colnames(oc) # dt.mnth Anomaly anomaly Solar variation
my_dat<-oc%>%dplyr::select(dt.mnth,anomaly)
library(Rssa)
L= 150
s_trend <- ssa(my_dat$anomaly, L = 150)
df_values.trd <- tibble(
  Component = seq_along(s_trend$sigma),
  SingularValue = s_trend$sigma
)

print(
  ggplot(df_values.trd[1:20, ], aes(Component, SingularValue)) +
    geom_point(color = "steelblue", size = 2) +
    geom_line(color = "steelblue") +
    scale_y_log10() +
    scale_x_continuous(breaks = 1:20) +
    labs(title = "SSA trend singular values (log scale)", x = "Component index", y = "Sigma") +
    theme_minimal(base_size = 11)
)
w_mat.trd <- as.matrix(wcor(s_trend, groups = 1:10))
print(
  expand.grid(F1 = 1:10, F2 = 1:10) %>%
    mutate(value = abs(as.vector(w_mat.trd))) %>%
    ggplot(aes(F1, F2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
    scale_x_continuous(breaks = 1:10, expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:10, expand = c(0, 0)) +
    coord_fixed() +
    labs(title = "W-trend correlation matrix", x = "Eigenvector index", y = "Eigenvector index",
         fill = "w-cor") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank())
)
rec.trend <- reconstruct(s_trend, groups = list(
  Global_Trend = c(1),
  natural= c(2:20)
))
print(
  tibble(dt.mnth = my_dat$dt.mnth,
         anomaly = my_dat$anomaly,
         Global = rec.trend$Global_Trend,
         natural = rec.trend$natural) %>%
    ggplot(aes(x = dt.mnth)) +
    geom_line(aes(y = anomaly, colour = "anomaly")) +
    geom_line(aes(y = Global, colour = "global"), linewidth = 1.2) +
    geom_line(aes(y = natural, colour = "long periods")) +
    labs(x = NULL, y = "anomaly (°C)", colour = NULL, title = "SSA-decomposed  trend") +
    theme_minimal(base_size = 11)
)
library(mgcv)
ssa_result=tibble(dt.mnth = my_dat$dt.mnth,
       anomaly = my_dat$anomaly,
       Global = rec.trend$Global_Trend,
       natural = rec.trend$natural)
ssa_result%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Global))+
  labs(x=NULL,title = "Global Trend ",subtitle = "ssa extracted")
ssa_result=ssa_result%>%
  mutate(nat.smth=predict(gam(natural ~ s(dt.mnth, k=17))))

ssa_result%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=natural,col="long.priods"))+
  geom_line(aes(y=nat.smth,col= "smth.priods"))+
  geom_line(aes(y=nat.smth),linetype = 2)+
  labs(x=NULL,title="Residuals of SSA extracted Trend")
