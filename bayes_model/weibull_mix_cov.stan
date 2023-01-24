data {
  // dimensions
  int<lower=0> N;             // number of observations
  int<lower=1> D;             // number of competing risks 
  vector[N] y;                // time for observation n
  int<lower=0> event[N];      // event status (1:D -- event type, 0:censor) for obs n
  int K;                      // number of covariates used for proportions 
  int M;                      // placeholder: the total number of used covariates [duplicates possible]
  matrix[N, K] Xprop;         // values for covariates used for proportions 
  int cov_dims[D,2];          // dimensions of covariates for each event type 
  matrix[N, M] X;             // covariate matrix [flattened] -- could use a ragged array once those are available 
}
transformed data {
  row_vector[K] zeros = rep_row_vector(0, K);
}

parameters {
  //simplex[D] theta;
  vector<lower=0>[D] alpha; 
  vector[D] mu;
  //vector[D] beta_raw[K];
  vector[M] beta;
  matrix[D-1, K] beta_pi_raw;
}

transformed parameters {
  matrix[D, K] beta_pi;
  matrix[N, D] lp;
  //vector[D] beta[K];
  beta_pi = append_row(beta_pi_raw, zeros);
  
  for (n in 1:N) {
    for (d in 1:D)
      lp[n,d] = mu[d] + dot_product(X[n,cov_dims[d,1]:cov_dims[d,2]], beta[cov_dims[d,1]:cov_dims[d,2]]);
  }
}

model {
  // priors
  vector[D] log_theta;
  alpha ~ gamma(1e-1, 1e-1);
  mu ~ normal(0, 5);
  to_vector(beta_pi_raw) ~ normal(0,1);
  to_vector(beta) ~ normal(0,1);

  // likelihood
  for (n in 1:N) {
      log_theta = log(softmax(beta_pi * transpose(Xprop[n,])));
      if (event[n] != 0){
          y[n] ~ weibull(alpha[event[n]], exp(-(lp[n,event[n]])/alpha[event[n]]));
          target += log_theta[event[n]];
      }
      else {
        vector[D] lps = rep_vector(0, D);
        for (d in 1:D){
          lps[d] = log_theta[d] + weibull_lccdf(y[n] | alpha[d], exp(-(lp[n,d])/alpha[d]));
        }
        target += log_sum_exp(lps);
      }
  }
}
generated quantities {
  matrix[N,D] pr = rep_matrix(0, N, D);
  
  for (n in 1:N) {
    vector[D] lps = log(softmax(beta_pi * transpose(Xprop[n,])));
    if (event[n]){
      pr[n, event[n]] = 1;
    }
    else{
      for (d in 1:D){
        lps[d] += weibull_lccdf(y[n] | alpha[d], exp(-(lp[n,d])/alpha[d]));
      }
    }
    pr[n,] = lps' - log_sum_exp(lps);
    pr[n,] = exp(pr[n,]);
  }
}
