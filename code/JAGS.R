# Header ------------------------------------------------------------
# Fitting a hierarchical linear model in JAGS
# Andrew Parnell and Ahmed Ali

# In this code we generate some data from a single level hierarchical model (equivalently a random effects model) and fit it using JAGS. We then interpret the output

# Some boiler plate code to clear the workspace, and load in required packages
rm(list = ls()) # Clear the workspace
library(R2jags)
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
mu = rnorm(n = 1, mu_mu, sd = sqrt(1/tau_mu))
tau = rgamma(n = 1, 1/alpha, beta)
omega = rgamma(n = 1, 1/alpha_omega, beta_omega)
mu_j = rnorm(n = M, mu, sd = sqrt(1/omega))

nj = sample(200:350, M, replace = TRUE) # Set the number of obs in each group between 5 and 10
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

# Jags code --------------------------------------------------------
# Jags code to fit the model to the simulated data

model_code = '
model
{
  # Likelihood
  for (i in 1:N) {
    y[i] ~ dnorm(mu_j[group[i]], tau^-1)
  }
  # Priors
  mu ~ dnorm(0, 1/0.5)
  omega ~ dgamma(0.2, 1)
  tau ~ dgamma(0.5, 1)
  for (j in 1:M) {
    mu_j[j] ~ dnorm(mu, omega^-1)
  }
}
'

# Set up the data
model_data = list(N = N, y = y, M = M, group = group)

# Choose the parameters to watch
model_parameters =  c("mu", "omega", "tau", "mu_j")

# Run the model
model_run = jags(data = model_data,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code))

model_run
# Simulated results -------------------------------------------------------

# Check the output - are the true values inside the 95% CI?
# Also look at the R-hat values - they need to be close to 1 if convergence has been achieved
plot(model_run)
print(model_run)
traceplot(model_run)

# Create a plot of the posterior mean regression line
post = print(model_run)
mu_j_pos = as.numeric(post$mean$mu_j) # Need the as.numeric otherwise it stores it as a weird 1D object
mu_j_pos = rep(mu_j_pos, nj)

data.frame(y = y, group = group, 
           muj = mu, mu_j_pos = mu_j_pos) %>% 
  ggplot(aes(y = y, x = group, group = group)) +
  geom_boxplot(fill = "#CD5C5C", alpha = 0.8) +
  geom_point(aes(y = mu), colour = '#0066ff', size = 3) +
  geom_point(aes(y = mu_j_pos), colour = '#ffff1a', size = 2) +
  ggtitle(label = "Blue points: true mean \nYellow points: posterior mean") +
  theme_bw(16)

# Blue (true) and red (predicted) points should be pretty close

# Real example ------------------------------------------------------------

# Load in
library(datasets)
head(women)

#Set up the data
jags_data = list(y = women[, 1],
                 N = dim(women)[1],
                 M = 2, 
                 group = rep(1:2, c(8, 7)))

# Plot
mus <- c(mean(jags_data$y[1:8]), mean(jags_data$y[9:15]))

# Set up jags model
jags_model=jags(jags_data,
                parameters.to.save = model_parameters,
                model.file = textConnection(model_code),
                n.chains = 4,
                n.iter=1000,
                n.burnin=200,
                n.thin=2)

# Plot output
print(jags_model)

# Create a plot of the posterior mean regression line
post = print(jags_model)
mu_j_pos = as.numeric(post$mean$mu_j) # Need the as.numeric otherwise it stores it as a weird 1D object

boxplot(jags_data$y ~ jags_data$group)
points(1:M, mus, col = 'blue', pch = 19)
points(1:M, mu_j_pos, col = 'red', pch = 19)
#--------------------------------------------------------------------

# Need better real data, better plots
