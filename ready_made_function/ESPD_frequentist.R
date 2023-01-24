# This script provides a set of custom functions that can be used to apply the Event-
# Specific Probabilities and Distributions (ESPD) approach for modelling two competing 
# events. This script provides the frequentist implementation of the ESPD approach. A
# Bayesian implementation is provided in another script. Although the functions can be  
# used for analyzing data for two competing events, they can be extended to account 
# for more than two competing events.
#
# The dcust(), pcust(), qcust(), and rcust() functions are custom functions for the 
# density, cumulative density, quantile, and sampling for a distribution. They are
# implemented so that these functions can be called for any of the below distributions
# by specifying the 'dist' argument, without the need to transform the estimated
# parameters to real scale. The code below estimates certain distribution parameters 
# on log-scale, in line with the 'flexsurv' package:
#   dist      p1                p2
# - gamma     shape = log       rate = log
# - gompertz  shape = real      rate = log
# - llogis    shape = log       scale = log
# - lnorm     meanlog = real    sdlog = log
# - weibull   shape = log       scale = log
dcust <- function(t, p1, p2, dist, logscale = FALSE) {
  d <- switch(
    dist,
    'gamma'     = dgamma    (x = t, shape = exp(p1), rate = exp(p2),  log = logscale),
    'gompertz'  = dgompertz (x = t, shape = p1,      rate = exp(p2),  log = logscale),
    'llogis'    = dllogis   (x = t, shape = exp(p1), scale = exp(p2), log = logscale),
    'lnorm'     = dlnorm    (x = t, meanlog = p1,    sdlog = exp(p2), log = logscale),
    'weibull'   = dweibull  (x = t, shape = exp(p1), scale = exp(p2), log = logscale),
    'weibullPH' = dweibullPH(x = t, shape = exp(p1), scale = exp(p2), log = logscale)
  )
  return(d)
}
pcust <- function(t, p1, p2, dist, lower.tail = TRUE, logscale = FALSE) {
  p <- switch(
    dist,
    'gamma'     = pgamma    (q = t, shape = exp(p1), rate = exp(p2),  lower.tail = lower.tail, log.p = logscale),
    'gompertz'  = pgompertz (q = t, shape = p1,      rate = exp(p2),  lower.tail = lower.tail, log.p = logscale),
    'llogis'    = pllogis   (q = t, shape = exp(p1), scale = exp(p2), lower.tail = lower.tail, log.p = logscale),
    'lnorm'     = plnorm    (q = t, meanlog = p1,    sdlog = exp(p2), lower.tail = lower.tail, log.p = logscale),
    'weibull'   = pweibull  (q = t, shape = exp(p1), scale = exp(p2), lower.tail = lower.tail, log.p = logscale),
    'weibullPH' = pweibullPH(q = t, shape = exp(p1), scale = exp(p2), lower.tail = lower.tail, log.p = logscale)
  )
  return(p)
}
qcust <- function(p, p1, p2, dist) {
  q <- switch(
    dist,
    'gamma'     = qgamma    (p = p, shape = exp(p1), rate = exp(p2)),
    'gompertz'  = qgompertz (p = p, shape = p1,      rate = exp(p2)),
    'llogis'    = qllogis   (p = p, shape = exp(p1), scale = exp(p2)),
    'lnorm'     = qlnorm    (p = p, meanlog = p1,    sdlog = exp(p2)),
    'weibull'   = qweibull  (p = p, shape = exp(p1), scale = exp(p2)),
    'weibullPH' = qweibullPH(p = p, shape = exp(p1), scale = exp(p2))
  )
  return(q)
}
rcust <- function(n, p1, p2, dist) {
  r <- switch(
    dist,
    'gamma'     = rgamma    (n = n, shape = exp(p1), rate = exp(p2)),
    'gompertz'  = rgompertz (n = n, shape = p1,      rate = exp(p2)),
    'llogis'    = rllogis   (n = n, shape = exp(p1), scale = exp(p2)),
    'lnorm'     = rlnorm    (n = n, meanlog = p1,    sdlog = exp(p2)),
    'weibull'   = rweibull  (n = n, shape = exp(p1), scale = exp(p2)),
    'weibullPH' = rweibullPH(n = n, shape = exp(p1), scale = exp(p2))
  )
  return(r)
}

# The getProb() and getLogit() functions are just convenience functions to
# transform between logits and probabilities.
getProb  <- function(logit) exp(logit) / (1 + exp(logit))
getLogit <- function(prob) log(prob / (1 - prob))

# The logLmix() function returns the log-likelihood for the mixture distribution
# of the competing events. This function is called by the optimization routine in 
# fitESPD() function to obtain the log-likelihood for a set of parameter values. 
# More details are provided within the function.
logLmix <- function(b, tt, ee, ww, dist1, dist2, X1, X11, X12, X21, X22) {
  
  # Input arguments:
  # - b       numeric vector of the coefficients that are being estimated
  # - tt      vector of times
  # - ee      vector of events
  # - ww      vector of weights
  # - dist1   character defining the distribution type for the first event
  # - dist2   character defining the distribution type for the second event
  # - X1      numeric covariate matrix: logit (p1) for experiencing event 1
  # - X11     numeric covariate matrix: distribution parameter 1 for event 1
  # - X12     numeric covariate matrix: distribution parameter 2 for event 1
  # - X21     numeric covariate matrix: distribution parameter 1 for event 2
  # - X22     numeric covariate matrix: distribution parameter 2 for event 2
  #
  # The following parameters define the mixture of competing events:
  # - p1    the probability of experiencing event 1
  # - p2    the probability of experiencing event 2 (1-p1)
  # - p11   the first parameter of the distribution for event 1
  # - p12   the second parameter of the distribution for event 1
  # - p21   the first parameter of the distribution for event 2
  # - p22   the second parameter of the distribution for event 2
  #
  # Note that what the distributions parameters are (e.g., shape, scale, rate,
  # meanlog, sdlog), and whether these are modeled on log- or real-scale, 
  # depends on the distribution type.
  
  # Because the number of coefficients may differ between use cases, they
  # are extracted from b based on the size of the covariate matrices. When
  # this function is called from the fitESPD() function, the order of p1, p11, 
  # p12, p21, p22 is always maintained.
  i_p1_start  <- 1
  i_p1_stop   <- ncol(X1)
  i_p11_start <- i_p1_stop + 1
  i_p11_stop  <- i_p11_start + ncol(X11) - 1
  i_p12_start <- i_p11_stop + 1
  i_p12_stop  <- i_p12_start + ncol(X12) - 1
  i_p21_start <- i_p12_stop + 1
  i_p21_stop  <- i_p21_start + ncol(X21) - 1
  i_p22_start <- i_p21_stop + 1
  i_p22_stop  <- i_p22_start + ncol(X22) - 1
  
  # Obtaining the individual-level distribution parameters
  p1  <- getProb(X1 %*% b[i_p1_start:i_p1_stop])
  p2  <- 1 - p1
  p11 <- X11 %*% b[i_p11_start:i_p11_stop]
  p12 <- X12 %*% b[i_p12_start:i_p12_stop]
  p21 <- X21 %*% b[i_p21_start:i_p21_stop]
  p22 <- X22 %*% b[i_p22_start:i_p22_stop]
  
  # Obtaining the (log-)likelihood using the custom distribution functions dcust()
  # and pcust() that call the appropriate distribution-specific functions, potentially
  # after transforming the parameters.
  L   <- (p1*dcust(tt, p11, p12, dist1))^(ee==1) * (p2*dcust(tt, p21, p22, dist2))^(ee==2) * (p1*pcust(tt, p11, p12, dist1, F) + p2*pcust(tt, p21, p22, dist2, F))^(ee==0)
  logL <- ww*log(L)
  
  return(sum(logL))
  
}

# The fitESPD() function is the main interface to all the other functions and
# can be called to start the procedure to estimate the event-specific probabilities
# and distributions (ESPD). More detail is provided in the function.
fitESPD <- function(t, e, w = NULL, i = NULL, events = c('progress', 'death'), dist1, dist2, X1 = NULL, X11 = NULL, X12 = NULL, X21 = NULL, X22 = NULL, nsim = 10^5, optmethod = 'BFGS', niter = 10^3) {
  
  # Input arguments:
  # - t       vector of times
  # - e       vector of events
  # - w       vector of weights
  # - i       logical vector to select the observations to be includes
  # - events  vector specifying the competing events in 'e' to determine
  #           which is event 1 and which is event 2, as well as what is the
  #           censoring event (i.e., the remaining third value in 'e')
  # - dist1   character defining the distribution type for the first event
  # - dist2   character defining the distribution type for the second event
  # - X1      optional numeric covariate matrix: logit (p1) for event 1
  # - X11     optional numeric covariate matrix: dist parameter 1 for event 1
  # - X12     optional numeric covariate matrix: dist parameter 2 for event 1
  # - X21     optional numeric covariate matrix: dist parameter 1 for event 2
  # - X22     optional numeric covariate matrix: dist parameter 2 for event 2
  # - nsim    optional numeric defining the number of samples to be used to 
  #           obtain a start value for the probability of event 1 (p1) based  
  #           on a cause-specific hazards model
  # - niter   optional numeric defining the maximum number of iterations for
  #           the optimization procedure
  
  if(!require(flexsurv)) stop("This function requires the 'flexsurv' package.")
  if(!require(maxLik))   stop("This function requires the 'maxLik' package.")
  
  # Setting default values if not provided
  n <- length(t)
  if(is.null(w)) w <- rep(x = 1, times = n)
  if(is.null(i)) i <- rep(x = TRUE, times = n)
  
  # If no covariate matrix is defined, a matrix is made with only an intercept
  # for each parameter
  if(is.null(X1))  X1  <- matrix(data = 1, nrow = n, ncol = 1, dimnames = list(NULL, 'Intercept'))
  if(is.null(X11)) X11 <- matrix(data = 1, nrow = n, ncol = 1, dimnames = list(NULL, 'Intercept'))
  if(is.null(X12)) X12 <- matrix(data = 1, nrow = n, ncol = 1, dimnames = list(NULL, 'Intercept'))
  if(is.null(X21)) X21 <- matrix(data = 1, nrow = n, ncol = 1, dimnames = list(NULL, 'Intercept'))
  if(is.null(X22)) X22 <- matrix(data = 1, nrow = n, ncol = 1, dimnames = list(NULL, 'Intercept'))
  
  p1_names <- c(gamma = 'logShape', gompertz = 'Shape', llogis = 'logShape', lnorm = 'Meanlog', weibull = 'logShape', weibullPH = 'logShape')
  p2_names <- c(gamma = 'logRate', gompertz = 'logRate', llogis = 'logScale', lnorm = 'logSdlog', weibull = 'logScale', weibullPH = 'logScale')
  
  colnames(X1)  <- paste0(events[1], '_logitTheta_', colnames(X1))
  colnames(X11) <- paste0(events[1], '_', p1_names[dist1], '_', colnames(X11))
  colnames(X12) <- paste0(events[1], '_', p2_names[dist1], '_', colnames(X12))
  colnames(X21) <- paste0(events[2], '_', p1_names[dist2], '_', colnames(X21))
  colnames(X22) <- paste0(events[2], '_', p2_names[dist2], '_', colnames(X22))
  
  t <- t[i]
  e <- e[i]
  w <- w[i]
  X1  <- X1[i, , drop = FALSE]
  X11 <- X11[i, , drop = FALSE]
  X12 <- X12[i, , drop = FALSE]
  X21 <- X21[i, , drop = FALSE]
  X22 <- X22[i, , drop = FALSE]
  n <- length(t)
  
  # Some basic checks
  if(!is.numeric(t) | any(t <= 0)) stop('Vector t needs to be a positive numeric')
  if(!is.numeric(w) | any(w <= 0)) stop('Vector w needs to be a positive numeric')
  if(!is.matrix(X1)  | !is.numeric(X1)  | nrow(X1) != n)  stop('Matrix X1 needs to be a numeric matrix with the number of rows equals to the length of vector t')
  if(!is.matrix(X11) | !is.numeric(X11) | nrow(X11) != n) stop('Matrix X11 needs to be a numeric matrix with the number of rows equals to the length of vector t')
  if(!is.matrix(X12) | !is.numeric(X12) | nrow(X12) != n) stop('Matrix X12 needs to be a numeric matrix with the number of rows equals to the length of vector t')
  if(!is.matrix(X21) | !is.numeric(X21) | nrow(X21) != n) stop('Matrix X21 needs to be a numeric matrix with the number of rows equals to the length of vector t')
  if(!is.matrix(X22) | !is.numeric(X22) | nrow(X22) != n) stop('Matrix X22 needs to be a numeric matrix with the number of rows equals to the length of vector t')
  if(length(e) != n | length(unique(e)) > 3) stop('Vector e needs to be of the same length as vector t and contain max. 3 unique values')
  if(!all(events %in% e)) stop('Competing events defined in argument events not all in vector e')
  
  # Transform the events to numbers, with 0 for censoring
  e[!(e %in% events)] <- '0'
  e[e == events[1]]   <- '1'
  e[e == events[2]]   <- '2'
  e <- as.integer(e)
  
  # Obtain start values for the intercepts of the distribution parameters using
  # a cause-specific hazards model using the 'flexsurv' package. Consequently, 
  # the scale of the parameters here are in line with that of the package.
  fit_cs_1 <- tryCatch(flexsurv::flexsurvreg(formula = Surv(t, e==1) ~ 1, dist = dist1), error = function(e) NULL)
  fit_cs_2 <- tryCatch(flexsurv::flexsurvreg(formula = Surv(t, e==2) ~ 1, dist = dist2), error = function(e) NULL)
  
  if(is.null(fit_cs_1)) stop('Estimation of the start values for the', dist1, 'distribution was unsuccessful.')
  if(is.null(fit_cs_2)) stop('Estimation of the start values for the', dist2, 'distribution was unsuccessful.')
  
  # Obtaining a start value for the intercept of the probability for event 1  
  sim_cs <- cbind(
    rcust(nsim, coef(fit_cs_1)[1], coef(fit_cs_1)[2], dist1),
    rcust(nsim, coef(fit_cs_2)[1], coef(fit_cs_2)[2], dist2)
  )
  
  # Define the vector of start values for the coefficients, with naive values
  # for the non-intercept coefficients and the probability on logit scale
  b_start <- setNames(object = c(
    getLogit(mean(sim_cs[,1] < sim_cs[,2])), rep(0, ncol(X1)-1),
    coef(fit_cs_1)[1], rep(0, ncol(X11)-1),
    coef(fit_cs_1)[2], rep(0, ncol(X12)-1),
    coef(fit_cs_2)[1], rep(0, ncol(X21)-1),
    coef(fit_cs_2)[2], rep(0, ncol(X22)-1)
  ), nm = c(colnames(X1), colnames(X11), colnames(X12), colnames(X21), colnames(X22)))
  
  # Maximum likelihood estimation
  fit_mix <- maxLik::maxLik(
    logLik = logLmix, start = b_start, method = optmethod, control = list(iterlim=niter),
    tt = t, ee = e, ww = w, dist1 = dist1, dist2 = dist2, 
    X1 = X1, X11 = X11, X12 = X12, X21 = X21, X22 = X22
  )
  
  # Output object
  ls_out <- list(
    convergence = fit_mix$code,
    iterations  = fit_mix$iterations,
    maximum     = fit_mix$maximum,
    estimates   = fit_mix$estimate,
    vcov        = tryCatch(vcov(fit_mix), error = function(e) NULL),
    AIC         = AIC(fit_mix),
    start       = b_start,
    maxLik      = fit_mix
  )
  
  # Only return a result if the optimization converged
  if(ls_out$convergence == 0) {
    cat('\n');
    cat('Parameters of the', dist1, '-', dist2, 'mixture were successfully estimated in', ls_out$iterations, 'iterations:'); cat('\n'); 
    print(ls_out$estimates); cat('\n');
    cat('AIC:', ls_out$AIC); cat('\n'); cat('\n');
    return(ls_out)
  } else {
    cat('\n'); cat('Parameters of the', dist1, '-', dist2, 'mixture were unsuccessfully estimated.'); cat('\n'); 
  }
  
}
