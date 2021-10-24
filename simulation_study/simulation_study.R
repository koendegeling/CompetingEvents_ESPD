# This script demonstrates the application of the frequentist implementation of the 
# "Event-Specific Probabilities and Distributions" (ESPD) approach for modelling
# two competing events based on covariates.
#
# The simulation considers a hypothetical scenario of two competing events:
# progression and death. These two events will be modelled based on age and disease
# stage as covariates. The probability of the event being progression and the scale
# parameter for the Weibull time-to-progression distribution will be modelled based 
# on the disease stage, whereas the rate parameter for the Gompertz time-to-death
# distribution will be modelled based on age. The simulation is implemented in a 
# function so that it can be run for different sample sizes and censoring rates.
#
# The below code was written and tested for the following versions:
# - R           v4.0.3
# - data.table  v1.14.0
# - flexsurv    v2.0


### 1. INITIALIZATION ----

# Uncomment to clear the Global Environment and console
#rm(list = ls()); gc(); cat('\14')

# Packages
library(data.table)   # efficient data wrangling
library(flexsurv)     # survival analysis functions

# Custom function for the frequentist implementation of the ESPD approach
source(file = 'R_functions/ESPD_frequentist.R')



### 2. DATA SIMULATION FUNCTION ----

simESPD <- function(pars, n_sim = 500, censoring_rate = NULL, seed = 123) {
  
  set.seed(seed)
  
  # Simulate the population with the covariates
  df <- data.table(ID = 1:n_sim)[ , `:=` (
    age   = rnorm(n = .N, mean = 60, sd = 5) - 60,
    stage = sample(x = c('IA', 'IB', 'II'), size = .N, replace = TRUE)
  )]
  
  # Create the covariance matrices
  X_stage <- cbind('Intercept' = rep(1, nrow(df)), 'IB' = as.integer(df$stage == 'IB'), 'II' = as.integer(df$stage == 'II'))
  X_age   <- cbind('Intercept' = rep(1, nrow(df)), 'Age' = df$age)
  
  # Calculate the individual-level distribution parameters
  df[ , `:=` (
    prob1  = getProb(X_stage %*% pars$progress_logitTheta),
    shape1 = exp(pars$progress_logShape),
    scale1 = exp(X_stage %*% pars$progress_logScale),
    shape2 = pars$death_Shape,
    rate2  = exp(X_age %*% pars$death_logRate)
  )]
  
  # Sampling the event and time-to-events
  df[ , `:=` (
    event = fifelse(runif(.N) < prob1, 'progress', 'death'),
    time1 = rweibull(.N, shape1, scale1),
    time2 = rgompertz(.N, shape2, rate2)
  )]
  df[ , time := fifelse(event=='progress', time1, time2)]
  
  # Sample a censoring time, if censor_rate provided
  if(!is.null(censoring_rate)) {
    df[ , time_cens := rexp(.N, rate = censoring_rate)]
    df[ , `:=` (
      event = fifelse(time_cens < time, 'censored', event),
      time  = fifelse(time_cens < time, time_cens, time)
    )]
    cat('\nProportion of observations that are censored:', mean(df$event == 'censored'), '\n\n')
  }
  
  out <- list(
    time    = df$time,
    event   = df$event,
    X_age   = X_age,
    X_stage = X_stage
  )
  
  return(out)
  
}




### 3. SIMULATIONS ----

# True parameter values
pars <- list(
  progress_logitTheta = c('Intercept' = -0.4, 'IB' = 0.4, 'II' = 0.8),
  progress_logShape   = c('Intercept' = 0.7),
  progress_logScale   = c('Intercept' = 2, 'IB' = -0.2, 'II' = -0.6),
  death_Shape         = c('Intercept' = 0.1),
  death_logRate       = c('Intercept' = -3.5, 'Age' = 0.1)
)

# Simulate datasets
df_0  <- simESPD(pars = pars, n_sim = 10^4)
df_10 <- simESPD(pars = pars, n_sim = 10^4, censoring_rate = 0.0135)
df_30 <- simESPD(pars = pars, n_sim = 10^4, censoring_rate = 0.0459)
df_60 <- simESPD(pars = pars, n_sim = 10^4, censoring_rate = 0.1472)

# Naive illustration
fitESPD(t = df_0$time, e = df_0$event, events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')

# Fitting ESPD parameters
fit_0  <- fitESPD(t = df_0$time,  e = df_0$event,  X1 = df_0$X_stage,  X12 = df_0$X_stage,  X22 = df_0$X_age,  events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')
fit_10 <- fitESPD(t = df_10$time, e = df_10$event, X1 = df_10$X_stage, X12 = df_10$X_stage, X22 = df_10$X_age, events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')
fit_30 <- fitESPD(t = df_30$time, e = df_30$event, X1 = df_30$X_stage, X12 = df_30$X_stage, X22 = df_30$X_age, events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')
fit_60 <- fitESPD(t = df_60$time, e = df_60$event, X1 = df_60$X_stage, X12 = df_60$X_stage, X22 = df_60$X_age, events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')

# Comparison of coefficients
rbind(
  truth  = c(pars$progress_logitTheta, pars$progress_logShape, pars$progress_logScale, pars$death_Shape, pars$death_logRate),
  fit_0  = fit_0$estimates, 
  fit_10 = fit_10$estimates, 
  fit_30 = fit_30$estimates, 
  fit_60 = fit_60$estimates 
)

# Probability of progression
rbind(
  truth  = getProb(cumsum(pars$progress_logitTheta)),
  fit_0  = getProb(cumsum(fit_0$estimates[1:3])),
  fit_10 = getProb(cumsum(fit_10$estimates[1:3])),
  fit_30 = getProb(cumsum(fit_30$estimates[1:3])),
  fit_60 = getProb(cumsum(fit_60$estimates[1:3]))
)

# Shape parameter for progression
rbind(
  truth  = exp(pars$progress_logShape),
  fit_0  = exp(fit_0$estimates[4]),
  fit_10 = exp(fit_10$estimates[4]),
  fit_30 = exp(fit_30$estimates[4]),
  fit_60 = exp(fit_60$estimates[4])
)

# Scale parameter for progression
rbind(
  truth  = exp(cumsum(pars$progress_logScale)),
  fit_0  = exp(cumsum(fit_0$estimates[5:7])),
  fit_10 = exp(cumsum(fit_10$estimates[5:7])),
  fit_30 = exp(cumsum(fit_30$estimates[5:7])),
  fit_60 = exp(cumsum(fit_60$estimates[5:7]))
)

# Shape parameter for death
rbind(
  truth  = pars$death_Shape,
  fit_0  = fit_0$estimates[8],
  fit_10 = fit_10$estimates[8],
  fit_30 = fit_30$estimates[8],
  fit_60 = fit_60$estimates[8]
)

# Rate parameter for death (40y)
rbind(
  truth  = exp(sum(c(1, 40) * pars$death_logRate)),
  fit_0  = exp(sum(c(1, 40) * fit_0$estimates[9:10])),
  fit_10 = exp(sum(c(1, 40) * fit_10$estimates[9:10])),
  fit_30 = exp(sum(c(1, 40) * fit_30$estimates[9:10])),
  fit_60 = exp(sum(c(1, 40) * fit_60$estimates[9:10]))
)

# Rate parameter for death (60y)
rbind(
  truth  = exp(sum(c(1, 60) * pars$death_logRate)),
  fit_0  = exp(sum(c(1, 60) * fit_0$estimates[9:10])),
  fit_10 = exp(sum(c(1, 60) * fit_10$estimates[9:10])),
  fit_30 = exp(sum(c(1, 60) * fit_30$estimates[9:10])),
  fit_60 = exp(sum(c(1, 60) * fit_60$estimates[9:10]))
)

# Rate parameter for death (80y)
rbind(
  truth  = exp(sum(c(1, 80) * pars$death_logRate)),
  fit_0  = exp(sum(c(1, 80) * fit_0$estimates[9:10])),
  fit_10 = exp(sum(c(1, 80) * fit_10$estimates[9:10])),
  fit_30 = exp(sum(c(1, 80) * fit_30$estimates[9:10])),
  fit_60 = exp(sum(c(1, 80) * fit_60$estimates[9:10]))
)



