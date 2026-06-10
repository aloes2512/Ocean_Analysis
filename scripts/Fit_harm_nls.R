library(tidyverse)
NOAA.source <- readRDS("NOAA.psl.rds")
library(zoo)
NOAA.psl <- NOAA.source$data%>%
  map(~ .x %>% dplyr::mutate(yr.mnth = as.yearmon(date))) %>%
  map(~ .x %>% dplyr::select(index, date, SST))


library(itsmr)
# select WPWP as example
mydat=NOAA.psl$WPWP%>%mutate(SST=smooth.fft(SST,0.01))
M4=c("season",12,"trend",4)
M0=c("season",12)
mydat.smth4=mydat%>%dplyr::mutate(yr.mnth = as.yearmon(date),
                                 res=Resid(SST,M4),
                                 trd=SST-res)
# season 12 removed=> fitted polynomial degree
my.poly4=mydat.smth4$trd
my.dff4=diff(my.poly4)
my.dff=my.dff4
#==================================
# example
ply=mydat.smth4$trd

# function to fit harmonic to polynomial
fitharm= function(ply){
  y=ply
  require(zoo)
  n=length(y)
  y.df=c(NA,diff(y))%>%na.locf()
  f.har=1/(4*(which.max(y.df)-which.min(y)))
  t0=which.max(y.df)-2*which.min(y) #226
  A= y[which.min(y)] # -0.74146
  y.fit=A*sin(2*pi*f.har*((1:n)-1+t0))
  return(y.fit)
}
# exa
tibble(t=mydat$date,y=fitharm(ply),y.pol=ply)%>%
  ggplot(aes(x=t,y=y))+
  geom_line(col=2)+
  geom_line(aes(y=y.pol))



t=as.yearmon(mydat$date)
tnm=as.numeric(t)
tibble(t = t,  my.dff = c(NA, my.dff4)) %>%
  ggplot(aes(x = t, y = my.dff)) +
  geom_line() +
  annotate("point", x = t[which.max(my.dff4)],
           y = my.dff[which.max(my.dff4)],
           colour = 2, size = 3) +
  annotate("point", x = t[which.min(my.poly4)],
           y = my.dff[which.min(my.poly4)],
           colour = 2, size = 3) +
  labs(title = "First Derivative of Poly4")
# adding  minimum poly and turning-point and text
tibble(t = t, tnm = tnm, my.dff = c(NA, my.dff4)) %>%
  ggplot(aes(x = t, y = my.dff)) +
  geom_line() +
  # Turning point
  annotate("point", x = t[which.max(my.dff4) + 1],
           y = my.dff4[which.max(my.dff4)],
           colour = 2, size = 3) +
  annotate("text", x = t[which.max(my.dff4) + 1] + diff(range(t))*0.02,
           y = my.dff4[which.max(my.dff4)] + diff(range(my.dff4))*0.05,
           label = "Turning Point Poly.4", colour = 2, size = 3) +
  # Minimum point
  annotate("point", x = t[which.min(my.poly4)],
           y = my.dff4[which.min(my.poly4)],
           colour = 2, size = 3) +
  annotate("text", x = t[which.min(my.poly4)] + diff(range(t))*0.02,
           y = my.dff4[which.min(my.poly4)] + diff(range(my.dff4))*0.05,
           label = "Minimum Poly.4", colour = 2, size = 3) +
  labs(title = "First Derivative of Poly4")

# estimate harmonic parms fitting the polynomial
quarter.prd=which.max(c(NA,my.dff4))-which.min(my.poly4)
prd=quarter.prd*4 # 4172 month  347 years 8 month (2/3 years)
## max gradient
my.dff4=diff(my.poly4)
t[which.max(my.dff4)] # "Oct 2008; 1858
length(t) # 2064
which.max(my.dff4)-2*quarter.prd # -228 month; -19 years
# sin function with parms A =amplitude; f.poly4= frequency; t0= shift of sin(0)
t0=(which.max(my.dff4)-2*quarter.prd)/12 #19 years
A=my.poly4[which.min(my.poly4)]
f.poly4=12/prd # 0.0002396932
phase.4= 2*pi*f.poly4*-(t0+tnm[1])
y.fit= A*sin(2*pi*f.poly4*(tnm-t0-tnm[1]))
#  or
y.fit1=A*sin(2*pi*f.poly4*tnm+phase.4)
#======
tibble(t=t,
       tnm=as.numeric(t),
       y.fit1=A*sin(2*pi*f.poly4*tnm+phase.4),
       y=mydat.smth4$trd)%>%
  ggplot(aes(x=t,y=y))+geom_line()+
  geom_line(aes(y=y.fit1),col=3)
#=================================

which.min(my.poly4) #816
t[which.min(my.poly4)] # Dec 1921
which.min(-A*sin(2*pi*f.poly4*(tnm+t0)))
mydat.no.ssn=mydat%>%dplyr::mutate(yr.mnth = as.yearmon(date),
                      res=Resid(SST,M0))
mydat.no.ssn%>%ggplot(aes(x=yr.mnth,y=res))+geom_line()
sgnl=mydat.no.ssn$res
tnm=as.numeric(mydat.no.ssn$yr.mnth)
fit3 <- lm(sgnl ~ stats::poly(tnm, 3))   # cubic trend
AIC(fit3) # -4107.728
library(broom)
y.fit=as.numeric(fitted(fit3))
tnm=as.numeric(mydat.smth4$yr.mnth)
tibble(t=as.numeric(mydat.smth4$yr.mnth)-1854,y=y.fit)%>%
  ggplot(aes(x=t,y=y))+geom_line()+
  geom_line(data = mydat.no.ssn,aes(x=t-t[1],y=res),col=2)
#=====
x=mydat.no.ssn$res
t <- seq_along(x)
harm_fn <- function(t, a, f, phi) {
   +a* sin(2 * pi * f * (t-t0) + phi)
}

fit_harm <- nls(x ~ harm_fn(t, a, f, phi),
                start = list(a = 1,f = 14*f.poly4, phi = pi))
AIC(fit_harm) # 870.3193
