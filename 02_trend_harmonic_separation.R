# 02_trend_harmonic_separation.R
# 02.1 SPLITTING BY FREQUENCY (periods < 9 years)
#splitting observed anomalies into low frequencies (periods> 10yrs)
#and high frequencies (periods < 9 years) as first step
#does NOT provide an unique date to separate
#modern times from pre-industrial
source("scripts/Anomalies_splitted_low_high.frequencies.R")
print(preind.plot)
# NOTE: THE RESULT DEPENDS ON THE SELECTED TIME 1950
#TO SEPARATE PRE-INDUSTRIAL
# --------
# Assuming that the RESIDUALS OF LOGISTIC TREND are an anomaly-model for pre-industrial times
# Rolling CORRELATION of Low-pass filtered ANOMALY and RESIDUALS OF LOGISTIC TREND with
#  a time span of 11 years do show a break of correlation (vertical line)
print(plt.cor)
# RESULT:
# 1.  the limit  separating pre-industrial and modern times
##  is @ 1987.833.. or as calendar date
## or as calendar date Oct 1987
# 2. the Ocean Anomaly values in pre-industrial times
##  are the residuals of a generalized logistic function
# 3. the long periods of pre-industrial times may be isolated
##  from the residuals
#  note: the method depends on several assumptions:
#        - the anomalies where not significant influenced
#          by industrial emissions before 1938
#        - the trend is modeled with sufficient accuracy using
#          generalized logistic function
#        - the logistic function is fitted with nonlinear model (nlm)
#          using as starting value the turning point of
#          trend fitted with polynomial 5th degree
#        - there may be artifacts resulting from the
#          time limit of the monthly observations 2116 month
#          from 1850 to 2116
rm(list=setdiff(ls(),"global.ts"))
#===================
#====================
# 02.2SPLITTING BY GAM_extrapolation and LOGISTIC FITTING
# split anomaly into pre-industrial times < 1938 and full obs. anomalies
# assumption industrial emissions had only marginal effect in ocean SST anomaly
# fit to gam model with 180 basis functions ~ 12 month separated
# extrapolate pre-ind to full times observed 2026.25
# industrial cause =difference of observed anomaly and xtrapol pre-ind
# difference == industrial caused
# fit difference with polynomial 5th degree.print(plt.gam.filt)
source("scripts/Gam_filter_xtrapol.R")
print(plt.gam.filt)
# ind. trend fitted with generalized logistic (Richards Curve) locked to the turning point
# using max.slope = estimated turning point to fit nlm logistics curve
print(plt.richards.trd)
#-------------------

#identified as pre-industrial and modern industrial times 1938 : 2026
# the pre-ind anomalies are superimposed harmonics and extended to 2026
# the difference to Anomaly extracts  natural changing anomaly
print(plt.rich.res)
# the plot shows also some ranges of superimposed harmonics:
##  sum ALL harmonics with periods < 45 years (perds<45)
##  all harmonics with periods between 18 and 45 years (18<perds< 30yr)
##  and harmonics with periods > 18 (perds >18yr)
#------------
# low-pass filtering Anomaly, Extended pre-industrial, and the difference of both
print(low.plt.decomposed)
# oberve that the difference is ascending mor rapidly than the observed Anomaly
# that means observed anomaly is reduced by the pre-ind trend or
# the modern times INDUSTRIAL effect shown as the difference was until ~ 2000
# BIGGER than OBSERVED ASCEND of global ocean anomaly
print(plt.rich.res)
rm(list=setdiff(ls(),c("trd.data_GAM","TP")))

#----------------------
#======================
# 02.3_Solar_locked
## A. variation of solar irradiation calculated from SILSO sunspot counts
##    trend solar locked fitted with splines and lm

source("scripts/Solar_locked_anomaly.R")
#Compare: solar-locked-trend polynomial-degree-7-trend
print(plt_trd7.trend.solar)
# polynomial 7th degree looks like a smoothed solar locked trend
# note: trend solar seems to end in a downward bending
# this could be an indication of an upper limit of the trend
# decompose solar-locked-trend residuals
print(plt.decomp_solar.lckd)
# note: extracted solar locked trend suggests:
#a secular trend superimposed by long periods
# decomposition:secular trend and long periods
source("scripts/DecomposeSolar.locked.trend.R")
print(plt.sol.lckd_trend)
# solar locked trend decomposed as:
#superposition of:
  ##clean secular trend
  ##long periods
  ## unknown known sum of noise and more periodical functions
# first identify secular trend

# ## Using LOESS to strip the smooth secular rise and
# separate the periodical functions of solar-locked-trend (trd.solar)
print(plt.sol.periods)
## show secular rise and deviation of solar locked trend
print(plt.sol.trend)
# summarize the solar locked decomposition into
## solar locked trend (trd.solar)
## periodic part of solar trend
## residuals of solar locked trend
# Analyze the Solar Locked Trend
source("scripts/Preind_Rolling_corr.sep.R")
# The solar locked trend was filtered using splines
# interpolating the intervals given by solar irradiation periods.
print(plt.sol.lckd)
# Trend isolated from the periodical changes (wriggles)
# gave the secular trend
# physical trend is assumed to stay within limits smaller than ± Inf.
# formulated mathematical it should be fitted to a
# function like tanh(time) or mor generally by Richards Curve.
# The generalized logistic to be fitted
#    Rch_t= A0+K*(1+exp(-B*(time-M)))^(-1/nu))
# depends on 5 parameters:
   ## an offset A0
   ## K related the distance of the upper and lower limit
   ## B determines the maximum slope
   ## M (time origin) the symmetry of the lower and upper part
   ## nu the asymetrie of the upper and lower slope
# fitting all 5 parms with one step is only possible
# if the observed data are trend towards saturation
#  if the timeseries does not map such a trend the solutions
#  for the 5 parms become ambiguous because the size of the slope depends
# depends on the length of the time-series
# to extract the trend in modern times fitted to a
# generalized logistic curve (Richards-Curve) therefore is done in two steps:
# A: estimating some value for M = 1970 and fitting the observed data to 4 parms:
print(plt.1st.iteration)
# B: linear fitting the 1st logistic curve to the observed timeseries

#========================
# SUMMARY TRENDS
# COMPARE 4  anomaly trends obtained with different filter:

              # GAM (gamtrd180),
              # low-pass (lowP.filt),
              # polynom 7th degree (trd7)
              # solar locked trend (trend.solar)
source("scripts/Extract_trend.R")
print(plt.trd)
rm(list=ls())
#==============
