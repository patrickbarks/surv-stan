functions {
  real expo_surv(real t, real l) {
    return exp(-l*t);
  }
}

data {
  int<lower=0> N;
  int<lower=0> N_pred;
  real t[N];
  vector[N_pred] t_new;
}

parameters {
  real<lower=0> l;
}

model {
  l ~ cauchy(0, 30);
  
  for(i in 1:N) {
    t[i] ~ exponential(l);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N_pred] surv_pred;
  
  for (i in 1:N) {
    log_lik[i] <- exponential_log(t[i], l);
  }
  
  for(j in 1:N_pred) {
    surv_pred[j] <- expo_surv(t_new[j], l);
  }
}
