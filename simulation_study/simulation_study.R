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
# - flexmix     v2.3-17
# - parallel    v4.0.3


## 1. INITIALIZATION ----

# Uncomment to clear the Global Environment and console
rm(list = ls()); gc(); cat('\14')

# Packages
library(data.table)   # efficient data wrangling
library(flexsurv)     # survival analysis functions
library(flexmix)      # Kullback-Leibler divergence
library(parallel)     # running analyses in parallel

# Custom function for the frequentist implementation of the ESPD approach
source(file = 'R_functions/ESPD_frequentist.R')



## 2. DATA SIMULATION FUNCTION ----

simESPD <- function(pars, n_sim = 500, seed = 123) {
  
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

  out <- list(
    time    = df$time,
    event   = df$event,
    X_age   = X_age,
    X_stage = X_stage
  )
  
  return(out)
  
}

censorEPSD <- function(obj, p_censoring = 0, start_censoring_rate = NULL, step_censoring_rate = 0.001, seed = 123) {
  
  set.seed(seed)
  
  df <- data.table(
    time  = obj$time,
    event = obj$event
  )
  
  # Censor the sample by increasing the censoring rate until the desired level is reached
  n_steps        <- 0
  p_censored     <- 0
  censoring_rate <- if(is.null(start_censoring_rate)) {step_censoring_rate} else {start_censoring_rate}
  while(p_censored < p_censoring) {
    
    n_steps <- n_steps + 1
    censoring_rate <- censoring_rate + step_censoring_rate
    
    df_cens <- copy(df)
    df_cens[ , time_cens := rexp(.N, rate = censoring_rate)]
    df_cens[ , `:=` (
      event = fifelse(time_cens < time, 'censored', event),
      time  = fifelse(time_cens < time, time_cens, time)
    )]
    
    p_censored <- mean(df_cens$event == 'censored')
    
    if(n_steps == 1 & p_censored > p_censoring) censoring_rate <- 0
    
  }
  
  cat('\nProportion of observations', round(p_censored, digits = 2), 'censored in', n_steps, 'steps using a censoring rate of', round(censoring_rate, digits = 4), '\n\n')
  
  out <- list(
    time    = df_cens$time,
    event   = df_cens$event,
    X_age   = obj$X_age,
    X_stage = obj$X_stage
  )
  
  return(out)  
}


## 3. TEST SIMULATIONS ----

# True parameter values
pars <- list(
  progress_logitTheta = c('Intercept' = -0.4, 'IB' = 0.4, 'II' = 0.8),
  progress_logShape   = c('Intercept' = 0.7),
  progress_logScale   = c('Intercept' = 2, 'IB' = -0.2, 'II' = -0.6),
  death_Shape         = c('Intercept' = 0.1),
  death_logRate       = c('Intercept' = -3.5, 'Age' = 0.1)
)

# Simulate datasets
ESPD_0  <- simESPD(pars = pars, n_sim = 10^4)
ESPD_10 <- censorEPSD(obj = ESPD_0, p_censoring = 0.1, start_censoring_rate = 0.01) #censoring_rate = 0.0135)
ESPD_30 <- censorEPSD(obj = ESPD_0, p_censoring = 0.30, start_censoring_rate = 0.04) #censoring_rate = 0.0459)
ESPD_60 <- censorEPSD(obj = ESPD_0, p_censoring = 0.60, start_censoring_rate = 0.14) #censoring_rate = 0.1472)

# Naive illustration
fitESPD(t = ESPD_0$time, e = ESPD_0$event, events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')

# Fitting ESPD parameters
fit_0  <- fitESPD(t = ESPD_0$time,  e = ESPD_0$event,  X1 = ESPD_0$X_stage,  X12 = ESPD_0$X_stage,  X22 = ESPD_0$X_age,  events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')
fit_10 <- fitESPD(t = ESPD_10$time, e = ESPD_10$event, X1 = ESPD_10$X_stage, X12 = ESPD_10$X_stage, X22 = ESPD_10$X_age, events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')
fit_30 <- fitESPD(t = ESPD_30$time, e = ESPD_30$event, X1 = ESPD_30$X_stage, X12 = ESPD_30$X_stage, X22 = ESPD_30$X_age, events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')
fit_60 <- fitESPD(t = ESPD_60$time, e = ESPD_60$event, X1 = ESPD_60$X_stage, X12 = ESPD_60$X_stage, X22 = ESPD_60$X_age, events = c('progress', 'death'), dist1 = 'weibull', dist2 = 'gompertz')

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




## 4. SIMULATION STUDY ----

getOutcomes <- function(t, e, TTP_min, TTP_max, TTD_min, TTD_max) {
  
  out <- list(
    prob_TTP = mean(e == 'progress'),
    mean_TTP = mean(t[e == 'progress']),
    mean_TTD = mean(t[e == 'death']),
    dens_TTP = density(x = t[e == 'progress'], n = 1024, from = TTP_min, to = TTP_max),
    dens_TTD = density(x = t[e == 'death'], n = 1024, from = TTD_min, to = TTD_max)
  )
  
  return(out)
  
}
getErrors <- function(sim, pop) {
  
  sim <- unname(sim)
  pop <- unname(pop)
  
  error <- sim - pop
  
  out <- c(
    E   = error,
    AE  = abs(error),
    RE  = error / pop,
    ARE = abs(error) / pop
  )
  
  return(out)
  
}
custSum <- function(x, ndigits = 2) {
  
  if(all(is.na(x)) | all(x == 0)) {
    out <- '-'
  } else {
    est <- mean(x, na.rm = TRUE)
    lb  <- quantile(x, 0.025, names = FALSE, na.rm = TRUE)
    ub  <- quantile(x, 0.975, names = FALSE, na.rm = TRUE)
    out <- paste0(formatC(est, digits = ndigits, format = 'f'), ' (', formatC(lb, digits = ndigits, format = 'f'), '; ', formatC(ub, digits = ndigits, format = 'f'), ')')
  }
  
  return(out)
  
}

runSIM <- function(pars, p_censoring, n_sample, n_runs = 1000, n_sim = 10^5, n_cores = 12, seed = 123) {
  
  set.seed(seed)
  n_scenarios <- length(p_censoring) * length(n_sample)
  i_scenario <- 0
  
  ls_out <- list()
  
  ESPD_pop <- simESPD(pars = pars, n_sim = n_sim)
  TTP_min  <- TTD_min <- 0
  TTP_max  <- max(ESPD_pop$time[ESPD_pop$event == 'progress'])
  TTD_max  <- max(ESPD_pop$time[ESPD_pop$event == 'death'])
  out_pop  <- getOutcomes(t = ESPD_pop$time, e = ESPD_pop$event) 
  
  for(i_censoring in p_censoring) {
    for(i_sample in n_sample) {
      
      i_scenario <- i_scenario + 1
      print(Sys.time())
      cat('Starting scenario', i_scenario, 'out of', n_scenarios, '\n\n')
      
      cl <- makeForkCluster(nnodes = n_cores)
      
      out_run <- parSapply(cl, 1:n_runs, function(i_run) {
        
        set.seed(seed + i_run)
        
        out_run <- NULL
        while(is.null(out_run)) {
          
          out_run <- tryCatch(expr = {
            
            i <- sample(x = n_sim, size = i_sample, replace = TRUE)
            
            ESPD_sample <- list(
              time    = ESPD_pop$time[i],
              event   = ESPD_pop$event[i],
              X_age   = ESPD_pop$X_age[i, ],
              X_stage = ESPD_pop$X_stage[i, ]
            )
            
            if(i_censoring > 0) {
              
              ESPD_sample <- censorEPSD(
                obj = ESPD_sample,
                p_censoring = i_censoring, 
                start_censoring_rate = switch(as.character(i_censoring), '0.1' = 0.01, '0.3' = 0.04, '0.6' = 0.10, NULL)
              )
              
            }
            
            fit <- fitESPD(t = ESPD_sample$time, e = ESPD_sample$event, 
                           dist1 = 'weibull', X1 = ESPD_sample$X_stage, X12 = ESPD_sample$X_stage,
                           dist2 = 'gompertz', X22 = ESPD_sample$X_age)
            
            if(fit$convergence != 0) stop('Estimation of ESPD unsuccessful')
            
            fit_pars <- list(
              progress_logitTheta = fit$estimates[1:3],
              progress_logShape   = fit$estimates[4],
              progress_logScale   = fit$estimates[5:7], 
              death_Shape         = fit$estimates[8],
              death_logRate       = fit$estimates[9:10]
            )
            sim <- simESPD(pars = fit_pars, n_sim = n_sim)
            
            if(any(is.infinite(sim$time))) stop('Inf value sampled')
            
            out_sim <- getOutcomes(
              t       = sim$time, 
              e       = sim$event, 
              TTP_min = TTP_min,
              TTP_max = TTP_max,
              TTD_min = TTD_min,
              TTD_max = TTD_max
            )
            
            error_prob <- getErrors(out_sim$prob_TTP, out_pop$prob_TTP)
            error_TTP  <- getErrors(out_sim$mean_TTP, out_pop$mean_TTP)
            error_TTD  <- getErrors(out_sim$mean_TTD, out_pop$mean_TTD)
            KLD_TTP    <- KLdiv(cbind(out_sim$dens_TTP$y, out_pop$dens_TTP$y))[2, 1] # population as reference
            KLD_TTD    <- KLdiv(cbind(out_sim$dens_TTD$y, out_pop$dens_TTD$y))[2, 1]
            
            c(probTTP = error_prob, meanTTP = error_TTP, meanTTD = error_TTD, TTP.KLD = KLD_TTP, TTD.KLD = KLD_TTD)
            
          }, error = function(e) NULL)
          
        }
        
        out_run
        
      })
      
      stopCluster(cl)
      
      summary_run <- apply(out_run, 1, custSum)
      
      ls_out[[length(ls_out) + 1]] <- c(p_censoring = i_censoring, n_sample = i_sample, summary_run)
      
    }
  }
  
  out <- do.call(rbind, ls_out)
  
  return(out)
  
}


sim <- runSIM(pars = pars, p_censoring = c(0, 0.1, 0.3, 0.6), n_sample = c(50, 100, 200, 500), n_runs = 10000, n_sim = 10^6)

saveRDS(object = sim, file = 'simulation_study/sim 20211226.RDS')
write.csv(x = sim, file = 'simulation_study/sim 20211226.csv', row.names = FALSE)

sim[, grepl('p_|n_|.ARE', colnames(sim))]



