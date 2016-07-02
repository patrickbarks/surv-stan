// https://github.com/dwcoder/StanSurvivalBoilerplate/

functions {
  real logist_log(real t, real l, real g, real s) {
    return log(l) + g*t + (-s-1)/s * log(1 + s*l/g * (exp(g*t)-1));
  }
  
  real logist_surv(real t, real l, real g, real s) {
    return (1 + s*l/g * (exp(g*t) - 1))^(-1/s);
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
  real<lower=0> g;
  real<lower=0> s;
}

model {
  l ~ cauchy(0, 30);
  g ~ cauchy(0, 30);
  s ~ cauchy(0, 30);

  for(i in 1:N) {
    t[i] ~ logist(l, g, s);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N_pred] surv_pred;
  
  for (i in 1:N){
    log_lik[i] <- logist_log(t[i], l, g, s);
  }
  
  for(j in 1:N_pred) {
    surv_pred[j] <- logist_surv(t_new[j], l, g, s);
  }
}
