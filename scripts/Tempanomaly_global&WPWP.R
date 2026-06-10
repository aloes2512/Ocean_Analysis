#1 prepare data tempanomaly
library(ncdf4)
library(tidyverse)
require(terra) # extract is masked by tidyr
ERSST.nc=nc_open("~/Desktop/Klima_Energiewende/Daten/gistemp1200_GHCNv4_ERSSTv5.nc")
ref_date<-as.Date("1800-01-01")
#time  Size:1754
#long_name: time
#units: days since 1800-01-01 00:00:00
#bounds: time_bnds
# variable: tempanomaly
time=ncvar_get(ERSST.nc,"time")
dates     <- ref_date + as.integer(time) # 1880-01-15 to 2026-02-15
ersst= terra::rast("~/Desktop/Klima_Energiewende/Daten/gistemp1200_GHCNv4_ERSSTv5.nc")
tempanom=unlist(global(ersst,"mean",na.rm=T)) # ersst  variable: tempanomaly
length(tempanom) #1754
y= tempanom
Y.new=tibble(t=seq_along(y)-1,
         y=tempanom)
#1 prepare data
NOAA.psl=readRDS("~/projects/Global_Temp/data/NOAA.psl.monthly.rds")
M= c("season",12) # linear trend already removed from SST
NOAA.ts <- map(
  NOAA.psl[-13],
  ~ .x %>% mutate(sst = Resid(value, M) %>% smooth.ma(12)))%>%
  map(~mutate(.x,yr.mnth=yearmon(year+(month-1)/12))%>%
        dplyr::select(yr.mnth,sst))
# select NOAA.ts$WPWP
Y= NOAA.ts$WPWP
colnames(Y) #yr.mnth sst
y=Y$sst # range : -0.5871719  0.9571583
t=seq_along(y)-1 # range t: 0 2063
# for a gridded sequence (length=50)
# the best fitting harmonics amplitude frequency and phase are calculated and the "best" selected by RSS minimum
N <- length(y) # 172 years
f_grid <- seq(0.1/N, 1/N, length.out=50)  # low freqs, adjust range

N <- length(y) # 172 years
f_grid <- seq(0.1/N, 1/N, length.out=50)
# Closure: captures y/t/f, accepts only params
sse_func_f <- function(params) {
  amp <- params[1] # essential for following up with DEoptim
  phase <- params[2] # essential for following up with DEoptim
  sum((y - amp * cos(2 * pi * f * t + phase))^2, na.rm = TRUE)
}
lower <- c(0, -pi)
upper <- c(1.5 * max(abs(y)), pi)

fits <- purrr::map_dfr(f_grid, ~{
  f <- .x
  sse_func_f <- function(params) sum((y - params[1]*cos(2*pi*f*t + params[2]))^2)
  fit <- DEoptim(sse_func_f, lower, upper, control = list(itermax = 200,steptol=0.01,trace=FALSE))  # fewer iters for speed
  tibble(f = f, amp = fit$optim$bestmem[1], phase = fit$optim$bestmem[2], sse = fit$optim$bestval)
})
grid_reslts=tibble( f= fits$f,
                    a=fits$amp,
                    phs=fits$phase,
                    bestv=fits$sse)%>%arrange(bestv)

five_best=grid_reslts%>%head(5)
best_periods=1/five_best$f
range(head(grid_reslts$f)) # 0.0003777092 0.0004666983
range(head(grid_reslts$a)) # 0.30  0.37
range(head(grid_reslts$phs)) #0.8160991 1.4785136

