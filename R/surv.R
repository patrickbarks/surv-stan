
# libraries and source
library(ggplot2)
library(dplyr)
library(rstan)
source('R/waic-fn.R')


# rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# read data (lifespan and other life history data for 214 fronds of common duckweed)
dat_url <- 'http://datadryad.org/bitstream/handle/10255/dryad.71127/BarksAndLaird_RawData_Phase1.csv'
dat <- read.csv(url(dat_url), stringsAsFactors = F)


# prepare data for stan
t <- sample(dat$lifespan, 100)
N <- length(t)
N_pred <- 100
t_new <- seq(0, max(t), length = N_pred)


# get emprical survivorship
SurvFn <- function(t_fail) {
  t <- 0:max(t_fail)
  p_surv <- sapply(t, function (x) length(which(t_fail >= x)) / length(t_fail))
  return(data.frame(t, p_surv))
}

df_surv <- SurvFn(t)


# stan parameter inits
inits_expo <- function() { list(l = runif(1, 0, 1))}
inits_weib <- function() { list(l = runif(1, 0, 1), b = runif(1, 0, 1))}
inits_gomp <- function() { list(l = runif(1, 0, 1), g = runif(1, 0, 1))}
inits_logi <- function() { list(l = runif(1, 0, 1), g = runif(1, 0, 1), s = runif(1, 0, 1))}


# fit stan models
FitStan <- function(file, init) {
  fit <- stan(
    file = file,
    data = list(t = t, N = N, N_pred = N_pred, t_new = t_new),
    init = init,
    iter = 3000,
    thin = 2,
    chains = 2
  )
}

fit_expo <- FitStan('stan/surv-exponential.stan', inits_expo)
fit_weib <- FitStan('stan/surv-weibull.stan', inits_weib)
fit_gomp <- FitStan('stan/surv-gompertz.stan', inits_gomp)
fit_logi <- FitStan('stan/surv-logistic.stan', inits_logi)


# extract fitted survivorship values
fitted_expo <- rstan::extract(fit_expo, pars = "surv_pred")$surv_pred
fitted_weib <- rstan::extract(fit_weib, pars = "surv_pred")$surv_pred
fitted_gomp <- rstan::extract(fit_gomp, pars = "surv_pred")$surv_pred
fitted_logi <- rstan::extract(fit_logi, pars = "surv_pred")$surv_pred


# combine fitted values from each model into tidy df
GetFitted <- function(model, fitted, t) {
  low <- apply(fitted, 2, function(x) quantile(x, probs = 0.025))
  med <- apply(fitted, 2, function(x) quantile(x, probs = 0.5))
  upp <- apply(fitted, 2, function(x) quantile(x, probs = 0.975))
  
  return(cbind.data.frame(low, med, upp) %>% mutate(t = t, model = model))
}

fitted_full <- rbind.data.frame(
  GetFitted('Exponential', fitted_expo, t_new),
  GetFitted('Weibull', fitted_weib, t_new),
  GetFitted('Gompertz', fitted_gomp, t_new),
  GetFitted('Logistic', fitted_logi, t_new)
)


# arrange models by number of parameters
fitted_full$model <- factor(
  fitted_full$model,
  levels = c('Exponential', 'Weibull', 'Gompertz', 'Logistic')
)


# calculate waic using function from Vehtari and Gelman (2014)
# http://www.stat.columbia.edu/~gelman/research/unpublished/waic_stan.pdf
waic_expo <- waic(fit_expo)$total[1] %>% round(1)
waic_weib <- waic(fit_weib)$total[1] %>% round(1)
waic_gomp <- waic(fit_gomp)$total[1] %>% round(1)
waic_logi <- waic(fit_logi)$total[1] %>% round(1)

df_waic <- data.frame(
  model = c('Exponential', 'Weibull', 'Gompertz', 'Logistic'),
  waic = c(waic_expo, waic_weib, waic_gomp, waic_logi)
) %>% mutate(label = paste('WAIC =', waic))


# plot
p1 <- ggplot(fitted_full) +
  geom_ribbon(aes(x = t, ymin = low, ymax = upp), fill = '#bdc9e1', alpha = 0.7) +
  geom_line(aes(t, med), size = 0.5, col = 'darkblue') +
  geom_point(data = df_surv, aes(t, p_surv), shape = 1, size = 2) +
  geom_text(data = df_waic, aes(label = label, 0, 0.004), hjust = 0, vjust = 0) +
  scale_y_log10(breaks = c(0.01, 0.1, 1), labels = c('0.01', '0.1', '1')) +
  coord_cartesian(ylim = c(0.004, 1)) +
  facet_wrap(~ model, nrow = 2) +
  xlab('Time (days)') + ylab('Survivorship') +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        text = element_text(size = 13),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.4, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .4, 0, 0, unit = 'cm')))

ggsave('img/surv-plots.png', p1, height = 5, width = 6, units = 'in', dpi = 300)



### alternate method to get fit lines in R instead of Stan
# extract model parameters
l_expo <- rstan::extract(fit_expo, pars = 'l')$l
l_gomp <- rstan::extract(fit_gomp, pars = 'l')$l
g_gomp <- rstan::extract(fit_gomp, pars = 'g')$g
l_logi <- rstan::extract(fit_logi, pars = 'l')$l
g_logi <- rstan::extract(fit_logi, pars = 'g')$g
s_logi <- rstan::extract(fit_logi, pars = 's')$s
l_weib <- rstan::extract(fit_weib, pars = "l")$l
b_weib <- rstan::extract(fit_weib, pars = "b")$b

# survivorship functions
SurvExpo <- function(l, t) exp(-l*t)
SurvWeib <- function(l, b, t) exp(-(l*t)^b)
SurvGomp <- function(l, g, t) exp(-(l/g) * (exp(g*t)-1))
SurvLogi <- function(l, g, s, t) (1 + s*l/g * (exp(g*t) - 1))^(-1/s)

# get fitted survivorship values
fitted_expo <- mapply(SurvExpo, l_expo, MoreArgs = list(t = t_new)) %>% t()
fitted_weib <- mapply(SurvWeib, l_weib, b_weib, MoreArgs = list(t = t_new)) %>% t()
fitted_gomp <- mapply(SurvGomp, l_gomp, g_gomp, MoreArgs = list(t = t_new)) %>% t()
fitted_logi <- mapply(SurvLogi, l_logi, g_logi, s_logi, MoreArgs = list(t = t_new)) %>% t()
