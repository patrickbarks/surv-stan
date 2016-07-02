// https://github.com/dwcoder/StanSurvivalBoilerplate/

functions {
  real gompertz_log(real t, real l, real g) {
    return log(l) + g * t - (l/g) * (exp(g*t)-1);
  }
  
  real gompertz_surv(real t, real l, real g) {
    return exp(-(l/g) * (exp(g*t)-1));
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
}

model {
  l ~ cauchy(0, 30);
  g ~ cauchy(0, 30);

  for(i in 1:N) {
    t[i] ~ gompertz(l, g);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N_pred] surv_pred;
  
  for (i in 1:N){
    log_lik[i] <- gompertz_log(t[i], l, g);
  }
  
  for(j in 1:N_pred) {
    surv_pred[j] <- gompertz_surv(t_new[j], l, g);
  }
}
