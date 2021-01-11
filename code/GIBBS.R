# Header ------------------------------------------------------------
# GIBBS code for fitting a hierarchical linear model
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

posteriors <- function(pars, type){
  if(type == 'mu'){
    omega <- pars$omega
    mu_bar <- mean(pars$mu_j)
    m <- pars$m
    mu_mu <- pars$mu_mu
    tau_mu <- pars$tau_mu
    
    mean_post <- (omega*mu_bar*m + mu_mu * tau_mu)/(tau_mu + omega * m)
    sd_post <- sqrt( 1 / (tau_mu + omega * m))
    post_sample <- rnorm(n = 1, mean = mean_post, sd = sd_post)
  }
  
  if(type == 'omega'){
    alpha_o <- pars$alpha_o
    mu_j <- pars$mu_j
    mu <- pars$mu
    m <- pars$m
    beta_o <- pars$beta_o
    term_mu <- sum((mu_j - mu)^2)
    
    alpha_post <- (m/2 + alpha_o)
    beta_post <- (term_mu/2 + beta_o)
    post_sample <- rgamma(n = 1, shape = alpha_post, scale = 1/beta_post)
  }
  
  if(type == 'tau'){
    alpha <- pars$alpha
    mu_j <- pars$mu_j
    y <- pars$y
    n <- length(y)
    beta <- pars$beta
    term_mu <- sum((mu_j - mu)^2)
    
    term_mu <- c()
    for(j in unique(group)){
      y_j <- y[group == j]
      term_mu[j] <- sum((y_j - mu_j[j])^2)
    }
    
    term_mu_f <- sum(term_mu)
    
    alpha_post <- (n/2 + alpha)
    beta_post <- (term_mu_f/2 + beta)
    post_sample <- rgamma(n = 1, shape = alpha_post, scale = 1/beta_post)
  }
  
  
  if(type == 'mu_j'){
    omega <- pars$omega
    n_j <- table(group)
    mu <- pars$mu
    tau <- pars$tau
    y <- pars$y
    
    y_bar_j <- c()
    for(j in unique(group)){
      y_bar_j[j] <- mean(y[group == j])
    }
    
    post_sample <- c()
    for(j in unique(group)){
      mean_post <- (n_j[j] * y_bar_j[j] * tau + mu * omega)/(n_j[j] * tau + omega)
      sd_post <- sqrt( 1 / (n_j[j] * tau + omega))
      post_sample[j] <- rnorm(n = 1, mean = mean_post, sd = sd_post)
    }
  }

  return(post_sample)  
  
}

# Starting Values ---
n_iter = 10000
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


pars = list(omega = current_omega, mu = current_mu,
            mu_j = current_mu_j, tau = current_tau,
            m = n_distinct(group), 
            mu_mu = mu_mu, tau_mu = tau_mu, 
            alpha = alpha, y = y, beta = beta, n = length(y),
            group = group, alpha_o = alpha_omega, 
            beta_o = beta_omega)


for(i in 2:n_iter){
  
  # Sample, save and update 
  # mu
  new_mu <- posteriors(pars, "mu")
  mus[i] <- new_mu
  pars$mu <- new_mu

  # mu
  new_omega <- posteriors(pars, "omega")
  omegas[i] <- new_omega
  pars$omega <- new_omega
  
  # tau
  new_tau <- posteriors(pars, "tau")
  taus[i] <- new_tau
  pars$tau <- new_tau
  
  # mu_j
  new_mu_j <- posteriors(pars, "mu_j")
  mu_js[[i]] <- new_mu_j
  pars$mu_j <- new_mu_j
  
}

results <- data.frame(
  iter = 1:(n_iter), 
  mus = mus, omegas = omegas, taus = 1/taus, 
  mu_1 = map_dbl(mu_js, 1), 
  mu_2 = map_dbl(mu_js, 2), 
  mu_3 = map_dbl(mu_js, 3), 
  mu_4 = map_dbl(mu_js, 4)
)
# taus are still very small?
# tau = should be 0.82
# omegas = should be 1.2...
# is stan inverted?
# The normal proposal works??
# Omega and tau are reversed?
# target omega = 0.86

new.lab <- as_labeller(c(mus = "mu", 
                         omegas = "omega", 
                         taus = "tau"), label_parsed)

results %>% 
  slice(-1) %>% 
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


# --------------------------------------------------------------------    
