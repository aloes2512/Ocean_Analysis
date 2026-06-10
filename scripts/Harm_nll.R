frequ=mx.prms[1] #0.00748336
phase=mx.prms[2] # 3.019336
time=(1:N)-1
observed=NINO.res$res
SD=sd(observed)
mu=mean(observed)
harmonic_nll <- function(frequency=frequ, phase, amplitude ) {

  predicted <- amplitude * cos(2 * pi * frequ * time + phase)
  residuals <- observed- predicted
  nLL= -sum(dnorm(residuals, mean = mu, sd = SD, log = TRUE))  # Gaussian likelihood
  return(nLL)
}
# start values

# Define starting values and bounds
## max from gridded
start_vals <- list( phase =mx.prms[2],amplitude=0.5)
lower_bounds <- c( phase = mx.prms[2]*0.9,amplitude=0.2)  # Lower bounds
upper_bounds <- c( phase =mx.prms[2]*1.1,amplitude=2)   # Upper bounds
time=(1:N)-1
observed=NINO.res$res
# Fit the model with bounded optimization
fit <- mle2(
  minuslogl = harmonic_nll,
  start = start_vals,   # named list for optimizer
  method = "L-BFGS-B",  # Optimizer that supports bounds
  lower = lower_bounds,
  upper = upper_bounds,
  data = list(time = time, observed = NINO.res$res, sigma = 1))

rslt.opt=coef(fit)
Harm.mdl=tibble(t=time,
                y.mdl=rslt.opt[2]* cos(2 * pi * frequ * time + rslt.opt[1]))