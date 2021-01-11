# Header ------------------------------------------------------------
# Fitting a hierarchical linear model to the 'repeated' data 
# Bruna Wundervald

# Some boiler plate code to clear the workspace, and load in required packages
rm(list = ls()) # Clear the workspace
library(tidyverse)

# Maths --------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y_{ij} = response variable for observation i = 1,...,n_j 
# in group j = 1,..,M.
# N = total number of observation = sum_j(n_j)
# mu_j = mean for each group j
# mu = mean for mu_j 
# tau = inverse variance of y_ij
# k = multiplication factor for the variance of mu_j 

# Likelihood:
# y_{ij} ~ N(mu_j, tau^{-1})
# Priors
# mu_j ~ N(mu, k(tau^{-1}))
# mu ~ N(0, tau_mu) 
# tau ~ Gamma(alpha, beta)
# alpha = 0.5, beta = 1, alpha = 0.2, beta = 1

# Simulate data -----------------------------------------------------
alpha = 0.5; beta = 1; mu_mu = 0; tau_mu = 0.3; k = 0.6


# Set the seed so this is repeatable
set.seed(2021)
# Some R code to simulate data from the above model
M = 4 # Number of groups
tau = rgamma(n = 1, 1/alpha, beta)
mu_main = rnorm(n = 1, mu_mu, sd = sqrt(tau_mu))
mu_j = rnorm(n = M, mu_main, sd = k * (1/tau))

nj = sample(200:5500, M, replace = TRUE) # Set the number of obs in each group between 5 and 10
N = sum(nj)
group = rep(1:M, times = nj)
mu = rep(mu_j, nj)
y = rnorm(N, mean = mu, sd = sqrt(1/tau))

data.frame(y = y, group = group, 
           muj = mu) %>% 
  ggplot(aes(y = y, x = group, group = group)) +
  geom_boxplot(fill = "#CD5C5C", alpha = 0.8) +
  geom_point(aes(y = mu), colour = '#0066ff', size = 3) +
  theme_bw(18)
ggsave(file = "book2/img/boxplot.png")

data.frame(y = y, group = group, muj = mu) %>% 
  group_by(group) %>% 
  summarise(var = var(y), 
            mean = mean(y)) %>% 
  mutate(muj = mu_j)

# Useful functions --------------------------------------------------

# Posteriors --- 
posteriors <- function(params, tau_p = NULL, 
                       mu_j_p = NULL,
                       mu_p = NULL){
  N <- params$N
  k <- params$k
  beta <- params$beta
  alpha <- params$alpha
  group <- params$group
  mu_j <- params$mu_j
  tau <- params$tau
  tau_mu <- params$tau_mu
  mu_main <- params$mu_main
  y <- params$y
  m <-  length(mu_j)
  
  if(!is.null(tau_p)){
  alpha_tau <- (N + m)/2 + alpha

  term_mu <- c()
  term_mu_j <- c()
  for(j in unique(group)){
    y_j <- y[group == j]
    term_mu[j] <- sum((y_j - mu_j[j])^2)
    term_mu_j[j] <- (mu_j[j] - mu_main)^2
  }
  
  beta_tau <- sum(term_mu)/2 + beta + sum(term_mu_j)/(2 * k)
  post <- dgamma(tau_p, alpha_tau, beta_tau)
  }
  
  
  if(!is.null(mu_j_p)){
    mu_j_p <-  rep(mu_j_p, length(nj))
    mean_mu <- c()
    var_mu <- c()

    for(j in unique(group)){
      y_bar_j <- mean(y[group == j]) 
      mean_mu[j] <- tau * ((mu_main/k) +  y_bar_j * nj[j])/(nj[j] + 1/k)
      var_mu[j] <- ((nj[j] + 1/k))^(-1) * tau
    }
    post <- dnorm(mu_j_p, mean = mean_mu, sd = sqrt(var_mu))
  } 
  
  if(!is.null(mu_p)){
    mean_mu <- (tau/k) * mean(mu_j) * m / (tau_mu + (tau/k)*m)
    var_mu <- (tau_mu + (tau/k)*m)^(-1)
    post <- dnorm(mu_p, mean = mean_mu, sd = sqrt(var_mu))
  } 
  
  return(list(post = post))
  
}

params <- list(
  N = N, mu_j  = mu_j, y = y, 
  group = group, beta = beta, 
  alpha = alpha, tau = tau, 
  tau_mu = tau_mu, 
  mu_main = mu_main, k = k
)


df <- data.frame(
  taus = seq(1, 2.5, length.out = 1000), 
  mu = seq(-1.5, 1.5, length.out = 1000)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    density_tau = posteriors(tau_p = taus, params = params),
    density_mu = posteriors(mu_j_p = mu, params = params),
    density_mu_main = posteriors(mu_p = mu, params = params))

df %>% 
  ggplot(aes(x = taus, y = unlist(density_tau))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = tau, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = expression(tau), y = "Posterior density", 
       title = expression("Posterior distribution for "~tau))
ggsave(file = "book2/img/post_tau.png")

df$density_mu %>% 
  map_df(as_tibble) %>% 
  mutate(ind_mu = rep(1:length(nj), 1000), 
         mu = rep(df$mu, each = 4)
  ) %>% 
  group_by(ind_mu) %>% 
  filter(value == max(value)) %>% 
  ungroup() %>% 
  arrange(ind_mu) %>% 
  mutate(mu_j = mu_j)

df$density_mu %>% 
  map_df(as_tibble) %>% 
  mutate(ind_mu = rep(1:length(nj), 1000), 
         mu = rep(df$mu, each = 4)) %>% 
  ggplot(aes(x = mu, y = value)) +
  geom_line() +
  theme_light(18) +
  facet_wrap(~ind_mu, scales = 'free') +
  geom_vline(data = data.frame(mu_j = mu_j, 
                               ind_mu = 1:4),
             aes(xintercept = mu_j), linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = expression(mu[j]), y = "Posterior density",
       title = expression("Posterior distributions for "~mu[j]))
ggsave(file = "book2/img/post_mu_j.png")

df %>% 
  ggplot(aes(x = mu, y = unlist(density_mu_main))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = mu_main, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = expression(mu), y = "Posterior density",
       title = expression("Posterior distribution for "~mu))
ggsave(file = "book2/img/post_mu.png")
