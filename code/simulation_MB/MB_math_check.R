# Header ------------------------------------------------------------
# Fitting a hierarchical linear model to the 'repeated' data 
# Bruna Wundervald

# Some boiler plate code to clear the workspace, and load in required packages
rm(list = ls()) # Clear the workspace
library(tidyverse)

# This is based on the latest piece of maths 
# Maths --------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y_{ij} = response variable for observation i = 1,...,n_j 
# in group j = 1,..,M.
# N = total number of observation = sum_j(n_j)
# mu_j = mean for each group j
# mu = mean for mu_j 
# tau = inverse variance of y_ij
# k1 = multiplication factor for the variance of mu
# k2 = multiplication factor for the variance of mu_j 

# Likelihood:
# y_{ij} ~ N(mu_j, tau^{-1})
# Priors
# mu_j ~ N(mu, k_2(tau^{-1}))
# mu ~ N(0, k_2(tau^{-1})) 
# tau ~ Gamma(alpha, beta)

# Simulate data -----------------------------------------------------
alpha = 0.5; beta = 1; mu_mu = 0; k1 = 1.5 ; k2 = 0.5

# Set the seed so this is repeatable
set.seed(2023)
# Some R code to simulate data from the above model
M = 9 # Number of groups
tau = rgamma(n = 1, 1/alpha, beta)
tau_mu_j = k1 * (1/tau) 
tau_mu = k2 * (1/tau) 
mu_main = rnorm(n = 1, mu_mu, sd = sqrt(tau_mu))
mu_j = rnorm(n = M, mu_main, sd = sqrt(tau_mu_j))

nj = sample(50:750, M, replace = TRUE) # Set the number of obs in each group between 5 and 10
N = sum(nj)
group = rep(1:M, times = nj)
mu = rep(mu_j, nj)
#M1 = model.matrix(mu ~ factor(group) - 1)
y = rnorm(N, mean = mu, sd = sqrt(1/tau))

# var_y = (
#   k1 * M1 %*% t(M1) + diag(1, length(mu), length(mu))) * (1/tau)
# 
# #var_y <- diag(diag(var_y))
# 
# y <- MASS::mvrnorm(n = 1, mu = rep(mu_main, N), Sigma = var_y)

data.frame(y = y, group = group, 
           muj = mu) %>% 
  ggplot(aes(y = y, x = group, group = group)) +
  geom_boxplot(fill = "#CD5C5C", alpha = 0.8) +
  geom_point(aes(y = mu), colour = '#0066ff', size = 3) +
  theme_bw(18)
# ggsave(file = "book2/img/boxplot.png")

data.frame(y = y, group = group, muj = mu) %>% 
  group_by(group) %>% 
  summarise(var = var(y), 
            mean = mean(y)) %>% 
  mutate(muj = mu_j)

# Useful functions --------------------------------------------------
# Posteriors --- 

params <- list(
  N = N, mu_j  = mu_j, y = y, 
  group = group, beta = beta, 
  alpha = alpha, tau = tau, 
  tau_mu = tau_mu, 
  nj = nj,
  mu_main = mu_main, k1 = k1, k2 = k2
  # det_psi = det_psi,
  # in_psi = in_psi, 
  # psi = psi
)

posteriors <- function(params, tau_p = NULL, 
                       mu_j_p = NULL,
                       mu_p = NULL){
  N <- params$N
  k1 <- params$k1
  k2 <- params$k2
  nj <- params$nj
  beta <- params$beta
  alpha <- params$alpha
  group <- params$group
  mu_j <- params$mu_j
  tau <- params$tau
  mu_main <- params$mu_main
  y <- params$y
  m <-  length(mu_j)

  
  if(!is.null(tau_p)){
    #alpha_tau <- (N + m)/2 + alpha
    alpha_tau <- (N + m + 1)/2 + alpha
    term_mu <- c()
    term_mu_j <- c()
    
    for(j in unique(group)){
      y_j <- y[group == j]
      term_mu[j] <- sum((y_j - mu_j[j])^2)
      term_mu_j[j] <- (mu_j[j] - mu_main)^2
    }
    
    
    beta_tau <- (sum(term_mu)/2 + beta + 
                   sum(term_mu_j)/(2 * k1) + (mu_main^2)/(2 * k2))

    # 8192.516
    post <- dgamma(tau_p, alpha_tau, beta_tau)
  }
  
  
  if(!is.null(mu_j_p)){
    mu_j_p <-  rep(mu_j_p, length(nj))
    mean_mu <- c()
    var_mu <- c()
    
    for(j in unique(group)){
      y_bar_j <- mean(y[group == j]) 
      mean_mu[j] <- ((mu_main/k1) +  y_bar_j * nj[j])/(nj[j] + 1/k1)
      var_mu[j] <- (1/(nj[j] + 1/k1))* tau
    }
    post <- dnorm(mu_j_p, mean = mean_mu, sd = sqrt(var_mu))
  } 
  
  if(!is.null(mu_p)){
    mean_mu <- (tau/k1) * mean(mu_j) * m / (tau *(m/k1 + 1/k2))
    #mean_mu <- ((1/k) * mean(mu_j))/ (tau_mu + (tau/k)*m)
    var_mu <- (tau *(m/k1 + 1/k2))^(-1)
    post <- dnorm(mu_p, mean = mean_mu, sd = sqrt(var_mu))
  } 
  
  return(list(post = post))
  
}

# psi <- (k + diag(N))
# det_psi <- det(psi)
# in_psi <- solve(psi)


df <- data.frame(
  taus = seq(tau - 1, tau + 1, length.out = 1000), 
  mu = seq(min(mu_j) - 1, max(mu_j) + 1, length.out = 1000)) %>% 
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
# ggsave(file = "book2/img/post_tau.png")

df$density_mu %>%
  map_df(as_tibble) %>%
  mutate(ind_mu = rep(1:length(nj), times = 1000),
         mu = rep(df$mu, each = M)
  ) %>%
  ungroup() %>%
  group_by(ind_mu) %>%
  filter(value == max(value)) %>%
  ungroup() %>%
  arrange(ind_mu) %>%
  mutate(mu_j = mu_j)
#  write.table("book2/mus.txt")

mu_lims <- mu_j %>% 
  map(~seq(.x -  0.3 * sqrt(k2 * (1/tau)), 
           .x + 0.3 * sqrt(k2 * (1/tau)) , length.out = 1000)) %>% 
  unlist()

df$density_mu %>% 
  map_df(as_tibble) %>% 
  mutate(ind_mu = rep(1:length(nj), 1000), 
         mu = rep(df$mu, each = M)) %>% 
  arrange(ind_mu) %>% 
  mutate(mu_lims = mu_lims) %>% 
  group_by(ind_mu) %>% 
  mutate(mu_min =  min(mu_lims), mu_max = max(mu_lims)) %>% 
  ungroup() %>% 
  filter(mu >= mu_min & mu <= mu_max) %>% 
  ggplot(aes(x = mu, y = value)) +
  geom_line() +
  theme_light(18) +
  facet_wrap(~ind_mu, scales = 'free') +
  geom_vline(data = data.frame(mu_j = mu_j, 
                               ind_mu = 1:M),
             aes(xintercept = mu_j), linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size = 12)) +
  labs(x = expression(mu[j]), y = "Posterior density",
       title = expression("Posterior distributions for "~mu[j]))
# ggsave(file = "book2/img/post_mu_j.png")

df %>% 
  ggplot(aes(x = mu, y = unlist(density_mu_main))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = mu_main, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  xlim(-3, 3) +
  labs(x = expression(mu), y = "Posterior density",
       title = expression("Posterior distribution for "~mu))
# ggsave(file = "book2/img/post_mu.png")

# Now using a smaller k1 ---------------------- 
alpha = 0.5; beta = 1; mu_mu = 0; k1 = 0.5 ; k2 = 0.2

# Set the seed so this is repeatable
set.seed(2023)
# Some R code to simulate data from the above model
M = 9 # Number of groups
tau = rgamma(n = 1, 1/alpha, beta)
tau_mu_j = k1 * (1/tau) 
tau_mu = k2 * (1/tau) 
mu_main = rnorm(n = 1, mu_mu, sd = sqrt(tau_mu))
mu_j = rnorm(n = M, mu_main, sd = sqrt(tau_mu_j))

nj = sample(50:750, M, replace = TRUE) # Set the number of obs in each group between 5 and 10
N = sum(nj)
group = rep(1:M, times = nj)
mu = rep(mu_j, nj)
#M1 = model.matrix(mu ~ factor(group) - 1)
y = rnorm(N, mean = mu, sd = sqrt(1/tau))

data.frame(y = y, group = group, 
           muj = mu) %>% 
  ggplot(aes(y = y, x = group, group = group)) +
  geom_boxplot(fill = "#CD5C5C", alpha = 0.8) +
  geom_point(aes(y = mu), colour = '#0066ff', size = 3) +
  theme_bw(18)
#ggsave(file = "book2/img/boxplot_sm.png")

data.frame(y = y, group = group, muj = mu) %>% 
  group_by(group) %>% 
  summarise(var = var(y), 
            mean = mean(y)) %>% 
  mutate(muj = mu_j)

params <- list(
  N = N, mu_j  = mu_j, y = y, 
  group = group, beta = beta, 
  alpha = alpha, tau = tau, 
  nj = nj, 
  tau_mu = tau_mu, 
  mu_main = mu_main, k1 = k1, k2 = k2
)

df <- data.frame(
  taus = seq(tau - 1, tau + 1, length.out = 1000), 
  mu = seq(min(mu_j) - 1, max(mu_j) + 1, length.out = 1000)) %>% 
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
#ggsave(file = "book2/img/post_tau_ms.png")

df$density_mu %>% 
  map_df(as_tibble) %>%
  mutate(ind_mu = rep(1:length(nj), times = 1000), 
         mu = rep(df$mu, each = M)
  ) %>% 
  ungroup() %>% 
  group_by(ind_mu) %>% 
  filter(value == max(value)) %>% 
  ungroup() %>% 
  arrange(ind_mu) %>% 
  mutate(mu_j = mu_j) %>% 
  write.table("book2/mus_ms.txt")


mu_lims <- mu_j %>% 
  map(~seq(.x -  0.5 * sqrt(k1 * (1/tau)), 
           .x + 0.5 * sqrt(k1 * (1/tau)) , length.out = 1000)) %>% 
  unlist()

df$density_mu %>% 
  map_df(as_tibble) %>% 
  mutate(ind_mu = rep(1:length(nj), 1000), 
         mu = rep(df$mu, each = M)) %>% 
  arrange(ind_mu) %>% 
  mutate(mu_lims = mu_lims) %>% 
  group_by(ind_mu) %>% 
  mutate(mu_min =  min(mu_lims), mu_max = max(mu_lims)) %>% 
  ungroup() %>% 
  filter(mu >= mu_min & mu <= mu_max) %>% 
  ggplot(aes(x = mu, y = value)) +
  geom_line() +
  theme_light(18) +
  facet_wrap(~ind_mu, scales = 'free') +
  geom_vline(data = data.frame(mu_j = mu_j, 
                               ind_mu = 1:M),
             aes(xintercept = mu_j), linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size = 12)) +
  labs(x = expression(mu[j]), y = "Posterior density",
       title = expression("Posterior distributions for "~mu[j]))
#ggsave(file = "book2/img/post_mu_j_ms.png")

df %>% 
  ggplot(aes(x = mu, y = unlist(density_mu_main))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = mu_main, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  xlim(-2, 2) +
  labs(x = expression(mu), y = "Posterior density",
       title = expression("Posterior distribution for "~mu))
# ggsave(file = "book2/img/post_mu_ms.png")

# Plot data densities using posteriors -----
# map_mu_j <- df$density_mu %>% 
#   map_df(as_tibble) %>% 
#   mutate(ind_mu = rep(1:length(nj), 1000), 
#          mu = rep(df$mu, each = M)) %>% 
#   group_by(ind_mu) %>% 
#   filter(value == max(value)) %>% 
#   arrange(ind_mu)
# 
# map_tau <- data.frame(value = df$density_tau %>%
#   unlist(), 
#   tau = df$taus) %>% 
#   filter(value == max(value))
# 
# map_mu <- data.frame(value = df$density_mu %>%
#                         unlist(), 
#                       mu = df$mu) %>% 
#   filter(value == max(value))
#   
