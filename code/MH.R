# Header ------------------------------------------------------------
# MH code for fitting a hierarchical linear model
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
# omega = inverse variance of mu_j

# Likelihood:
# y_{ij} ~ N(mu_j, tau^{-1})
# Priors
# mu_j ~ N(mu, omega^{-1})
# mu ~ N(mu_mu, tau_mu) 
# tau ~ Gamma(alpha, beta)
# omega ~ Gamma(alpha_omega, beta_omega)
# alpha = 0.5, beta = 1, alpha_omega = 0.2, beta_omega = 1, 
# mu_mu = 0, tau_mu = 0.5 

# Simulate data -----------------------------------------------------
alpha = 0.5; beta = 1; alpha_omega = 0.2; beta_omega = 1; 
mu_mu = 0; tau_mu = 0.5

# Set the seed so this is repeatable
set.seed(123)
# Some R code to simulate data from the above model
M = 4 # Number of groups
mu_main = rnorm(n = 1, mu_mu, sd = sqrt(1/tau_mu))
tau = rgamma(n = 1, 1/alpha, beta)
omega = rgamma(n = 1, 1/alpha_omega, beta_omega)
mu_j = rnorm(n = M, mu_main, sd = sqrt(1/omega))

nj = sample(200:350, M, replace = TRUE) # Set the number of obs in each group between 5 and 10
N = sum(nj)
group = rep(1:M, times = nj)
mu = rep(mu_j, nj)
y = rnorm(N, mean = mu, sd = sqrt(1/tau))

# Useful functions --------------------------------------------------
likelihood <- function(y, group, pars){
  
  tau <- pars$tau
  mu_j <- pars$mu_j
  
  ds <- vector(mode = "numeric", length = n_distinct(group))
  for(j in unique(group)){
    y_j <- y[group == j]
    ds[j] <- sum(
      dnorm(y_j, mean = mu_j[j], sd = sqrt(1/tau), log = TRUE))
  }
  summll <-  sum(ds)
  return(summll)
}

priors <- function(pars, alpha_p = alpha, 
                   beta_p = beta, 
                   alpha_o_p = alpha_omega, 
                   beta_o_p = beta_omega, 
                   mu_mu_p = mu_mu, 
                   tau_mu_p = tau_mu){
  mu_j <- pars$mu_j
  tau <- pars$tau
  omega <- pars$omega
  mu <- pars$mu 
  
  tau_prior <-  dgamma(tau, 1/alpha_p, beta_p, log = TRUE)
  omega_prior <-  dgamma(omega, 1/alpha_p, beta_o_p, log = TRUE)
  mu_prior <-  dnorm(mu, mean = mu_mu_p, sd = sqrt(1/tau_mu_p), log = TRUE)
  mu_j_prior <-  dnorm(mu_j, mu, sd = sqrt(1/omega), log = TRUE)
  all_priors <- sum(tau_prior, omega_prior, 
                    mu_prior, mu_j_prior)
  return(all_priors)
}


posterior <- function(y, group, pars){
  mu_j <- pars$mu_j
  tau <- pars$tau
  omega <- pars$omega
  mu <- pars$mu 
  
  post <- 
    priors(pars = 
             list(mu_j = mu_j, tau = tau, omega = omega, mu = mu)) +
    likelihood(y, group, list(tau = tau, mu_j = mu_j))
  
  return(post)
  
}

q_theta <- function(old, new, alpha, beta = 0.2){
  old_q <- dgamma(old, shape = alpha/beta, scale = beta, log = TRUE)
  new_q <- dgamma(new, shape = alpha/beta, scale = beta, log = TRUE)
  
  return(list(old_q = old_q, new_q = new_q))
  
}


# Starting Values ---
n_iter = 2000
taus = mus = omegas = c()
mu_js = list()

current_mu_j = mu_j 
current_mu = mu_main
current_tau = tau
current_omega = omega

taus[1] = current_tau
mus[1] = current_mu
omegas[1] = current_omega
mu_js[[1]] = mu_j

acc_mu = 0
acc_tau = 0
acc_omega = 0
acc_mu_js = list(0, 0, 0, 0)


for(i in 1:n_iter){
  
  # Updating mu_j
  for (j in unique(group)) {
    
    mu_j_candidate = current_mu_j 
    mu_j_candidate[j] = rnorm(1, current_mu_j[j], sd = 0.5)
    
    #mu_all_new[j] = mu_j_new #change jth element
    new_post = posterior(
      y, group, 
      list(mu_j = mu_j_candidate, 
           tau = current_tau, mu = current_mu, omega = current_omega))
    
    old_post = posterior(
      y, group, 
      list(mu_j = current_mu_j, 
           tau = current_tau, mu = current_mu, omega = current_omega))
    
    u = runif(1)
    a = exp(new_post - old_post) 
    
    if(a > u){
      current_mu_j = mu_j_candidate
      acc_mu_js[[j]] = acc_mu_js[[j]] + 1
    } 
  }
  
  mu_js[[i+1]] = current_mu_j
  
  
  # Sampling tau - only change that and accept/reject
  #tau_candidate = abs(rnorm(1, current_tau, 0.2)) # A normal proposal??
  tau_candidate = rgamma(1, shape = current_tau/0.2, scale = 0.2)
  
  
  if(tau_candidate == 0){
    tau_candidate = tau_candidate + abs(rnorm(1, 0, 0.1))
  }
  
  qs_tau <- q_theta(old = current_tau, 
                    new = tau_candidate, alpha = current_tau)
  
  # qs_tau <- list(
  #   old_q = dgamma(current_tau, shape = tau_candidate/0.2, scale = 0.2, 
  #                  log = TRUE),
  #   new_q = dgamma(tau_candidate, shape = current_tau/0.2, scale = 0.2,
  #                  log = TRUE)
  #   
  # )
  
  
  new_post = posterior(
    y, group, 
    pars = list(tau = tau_candidate, mu_j = current_mu_j, 
                omega = current_omega, mu = current_mu))
  
  old_post = posterior(
    y, group, 
    list(mu_j = current_mu_j, tau = current_tau, 
         mu = current_mu, omega = current_omega))
  
  u  = runif(1)
  a = exp((new_post + qs_tau$old_q) - (old_post + qs_tau$new_q)) 
  
  if(a > u){
    current_tau = tau_candidate
    acc_tau = acc_tau +1
  } 
  
  taus <- c(taus, current_tau)
  
  
  # Sampling mu only 
  mu_candidate = rnorm(1, current_mu, sd = 0.5)
  
  new_post = posterior(
    y, group, 
    pars = list(tau = current_tau, mu_j = current_mu_j, 
                omega = current_omega, mu = mu_candidate))
  
  old_post = posterior(
    y, group, 
    list(mu_j = current_mu_j, tau = current_tau, 
         mu = current_mu, omega = current_omega))
  
  u  = runif(1)
  a = exp(new_post - old_post) 
  
  if(a > u){
    current_mu = mu_candidate
    acc_mu = acc_mu +1
  }
  
  mus <- c(mus, current_mu)
  
  # Sampling omega only 
  # omega_candidate = rgamma(1, shape = current_omega/2, scale = 2)
  # 
  # if(omega_candidate == 0){
  #  omega_candidate = omega_candidate + abs(rnorm(1, 0, 0.1))
  #  }
  
  #omega_candidate = abs(rnorm(1, current_omega, 0.2))
  
  omega_candidate = rgamma(1, shape = current_omega/0.2, scale = 0.2)
  
  qs_omega <- q_theta(current_omega, omega_candidate, alpha = current_omega)
  # qs_omega <- list(
  #   old_q = dgamma(current_omega, shape = omega_candidate/0.2, scale = 0.2, 
  #                  log = TRUE),
  #   new_q = dgamma(omega_candidate, shape = current_omega/0.2, scale = 0.2,
  #                  log = TRUE)
  #   
  # )
  
  
  new_post = posterior(
    y, group, 
    pars = list(tau = current_tau, mu_j = current_mu_j, 
                omega = omega_candidate, mu = current_mu))
  
  old_post = posterior(
    y, group, 
    list(mu_j = current_mu_j, tau = current_tau, mu = current_mu, 
         omega = current_omega))
  
  u  = runif(1)
  a = exp((new_post + qs_omega$old_q) - (old_post + qs_omega$new_q)) 
  
  if(a > u){
    current_omega = omega_candidate
    acc_omega = acc_omega +1
  }
  
  omegas <- c(omegas, current_omega)
}


c(acc_mu, acc_omega, acc_tau, acc_mu_js[[1]], 
  acc_mu_js[[2]], acc_mu_js[[3]], acc_mu_js[[4]])/n_iter

mu_j_candidate

results <- data.frame(
  iter = 1:(n_iter + 1), 
  mus = mus, omegas = omegas, taus = 1/taus, 
  mu_1 = map_dbl(mu_js, 1), 
  mu_2 = map_dbl(mu_js, 2), 
  mu_3 = map_dbl(mu_js, 3), 
  mu_4 = map_dbl(mu_js, 4)
)

new.lab <- as_labeller(c(mus = "mu", 
                         omegas = "omega", 
                         taus = "1/tau"), label_parsed)

results %>% 
  select(iter, mus, omegas, taus) %>% 
  gather(parameter, value, -iter) %>% 
  group_by(parameter) %>% 
  mutate(mean = mean(value)) %>% 
  ggplot(aes(x = iter, y = value), colour = 'grey80') + 
  geom_line() +
  geom_line(aes(y = mean), colour = "#CD5C5C", size = 1, 
            linetype = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  facet_wrap(~parameter, scales = 'free', nrow = 3, 
             labeller = new.lab) +
  ggtitle(label = 
            "Samples for each parameter") +
  xlab("Iteration") + ylab("Values") +
  theme_bw(16)

new.lab.mu <- as_labeller(c(mu_1 = "mu[1]", 
                            mu_2 = "mu[2]", 
                            mu_3 = "mu[3]", 
                            mu_4 = "mu[4]"), label_parsed)

results %>% 
  select(iter, mu_1, mu_2, mu_3, mu_4) %>% 
  gather(parameter, value, -iter) %>% 
  group_by(parameter) %>% 
  mutate(mean = mean(value)) %>% 
  ggplot(aes(x = iter, y = value), colour = 'grey80') + 
  geom_line() +
  geom_line(aes(y = mean), colour = "#CD5C5C", size = 1, 
            linetype = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  facet_wrap(~parameter, scales = 'free', nrow = 2, 
             labeller = new.lab.mu) +
  ggtitle(label = 
            expression("Samples for each  "~mu[j])) +
  xlab("Iteration") + ylab("Values") +
  theme_bw(16)

# --------------------------------------------------------------