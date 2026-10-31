# Rolling Cor
# 1 define reference non industrial periods as residuals
## logistic curve obtained from SSA_decomp_logistic
library(tidyverse)
library(minpack.lm)
Ocean_solar_anomaly<-readRDS("data/Ocean_solar_anomaly.rds")
dt.mnth=Ocean_solar_anomaly$dt.mnth


# calculate residuals using Rssa
library (Rssa)
my_signal=Ocean_solar_anomaly$anomaly
L <- 730
s <- ssa(my_signal, L = L)
# 1. Compute wcor
w <- wcor(s, groups = 1:20)
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
# 1. Extract the reconstructed trend component vector
rc <- reconstruct(s, groups = list(trend =1:2,
                                   Har2=4:5,
                                   Har3=6:9,
                                   Har4=10:20))
rc.full<- reconstruct(s,groups = list(all=1:40))
my_rconst<-tibble(dt.mnth=Ocean_solar_anomaly$dt.mnth,
       full=rc.full$all,
       trd=rc$trend,
       Har2=rc$Har2,
       Har3=rc$Har3)
my_rconst%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Har2))+
  geom_line(aes(y=Har3))


fit_log <- nlsLM(
  full ~ A / (1 + exp(-k * (dt.mnth - t0))) +b,
  data = my_rconst,
  start = list(A = 1.0, k = 0.03, t0 = 1990,b=0)
)
round(coef(fit_log),3) #  L       k    t0         b
                    #  1.41   0.046 2017.916 -0.128

x_infl=round(coef(fit_log)[[3]])

df_logstic=tibble(dt.mnth=dt.mnth,
       logstc=predict(fit_log))
P_infl=df_logstic%>%subset(dt.mnth==x_infl)

my_rconst=my_rconst%>%left_join(df_logstic)

my_rconst%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trd),col=3)+
  geom_line(aes(y=logstc),col=4)+
  geom_line(aes(y=full),col="grey")+
  geom_point(data = P_infl,aes(x=dt.mnth,y=logstc))+
  labs(x=NULL,y="anomaly",title="Logistic fitted to\nSSA-reconstructed anomaly ")
ggsave("figs/Logistic_SSA.anomaly.png")
# define difference of logistic and reconstructed ssa as pre.indres
## pre.indres = full-logstc
my_dat = my_rconst%>%
  mutate(pre.indres=full-logstc,vect_obs=full)
colnames(my_dat) # "dt.mnth""full""trd""Har2""Har3"
                 #"logstc""pre.indres" "vect_obs"

#==============
library(slider)
vector_obs<-my_dat$full
vector_res<-my_dat$pre.indres
rolling_corr2 <- slide_index_dbl(
  .x = 1:length(vector_obs),
  .i = 1:length(vector_obs),
  .f = function(idx) {
    # If the window doesn't have a full 11 years of data, return NA
    if(length(idx) < 132) return(NA)

    # Calculate simple vector correlation on the current slice of indices
    cor(vector_obs[idx], vector_res[idx])
  },
  .before = 132,   # Looks 131 steps backward + current step = 132 months
  .complete = TRUE # Automatically returns NA for incomplete windows at the start
)

RollingCor=tibble(dt.mnth=Ocean_solar_anomaly$dt.mnth,
                  rolling_corr2=unlist(rolling_corr2))#
# select break points
# plo shows cor==1 all times < 1958
cor.sel=RollingCor%>%
  subset( dt.mnth< 1958 )%>%
  pull(rolling_corr2)
# select cor min
break.cor=dt.mnth[which.min(cor.sel)] #1957.417
# the correlated values are 132 month 11 years
break1= break.cor-11 # 1946.417 start of breaking the cor is minus 11 years
cor.sel2=RollingCor%>%
  subset(dt.mnth< 1985 )%>%
  pull(rolling_corr2)
# total break of cor
break.<-dt.mnth[which.min(cor.sel2)] # 1977
break2<-RollingCor%>%subset(rolling_corr2>0.95)%>%
  arrange(desc(dt.mnth))%>%{.[1,]} # 1986
RollingCor%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=rolling_corr2))+
  geom_vline(xintercept = break1,linetype = 2)+
  geom_vline(xintercept = break2[[1]],linetype = 2)+
  labs(x=NULL,y="rolling correlation",
       title = "Rolling Correlation with 2 breaks",
       subtitle= "corr: 11 yrs interval of pre-industrial")

ggsave("figs/roll_cor_preind.png")
brk1<-break1
brk2<-round(break2[[1]])

my_dat%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=pre.indres,col= "pre-ind"))+
  geom_line(aes(y=full,col= "observ"),)+
  geom_line(aes(y=logstc,col="industr"))+
  geom_vline(
    xintercept = c(brk1, brk2),
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  )+
  geom_text(
    data = data.frame(x = c(brk1, brk2), label = c("brk1", "brk2")),
    aes(x = x, y = Inf, label = label), # y = Inf places them near the top of the plot area
    vjust = 1.5,                        # Adjust vertical position down from the top edge
    hjust = -0.2,                       # Offset slightly to the right of the line
    inherit.aes = FALSE,
    size = 4
  )+
  geom_point(data = P_infl,aes(x=dt.mnth,y=logstc))+
  labs(x=NULL,y="anomaly",
       title = "Industr. Effect on Global Anomaly",
       subtitle=paste0("transition period: ",brk1," to ",brk2))

ggsave("figs/Observed_Industr_separated.png")




