library(tidyverse)
library(reshape2)
library(mvtnorm)
library(gridExtra)
library(R2jags)
set.seed(2025)

# Model description -------------------------------------------------------
# We have R_{ij} ~ N(mu_j, tau^{-1})
# mu_j ~ N(mu, k1 tau^{-1}/P)
# mu ~ N(0, k2 tau^{-1}/P)
# tau ~ Ga(alpha, beta)
# --------------------------------------------------

# Number of obs
N <- n <- 1000 # You can increase this, but some of the comparisons below get very slow
J <- 10 # Number of groups

# Other hyper-parameters defined above
alpha <- 1
beta <- 1

# Simulate data -----------------------------------------------------------
# Create a factor with the groups
x <- factor(sample(1:J, size = n, replace = TRUE))
# Group sizes
n_j <- table(x)
# Group allocation matrix
# Mean shrinkage parameter (always fixed)
k2 <- 5
# Conditional simulation
k1_TRUE <- 8
tau_TRUE <- rgamma(1, alpha, beta)
y <- vector(length = n)
# Using 5 trees
P <- 5

mu_true <- vector(length = P)
mu_j_TRUE <- list()

# for(i in 1:P){
#   mu_TRUE[i] <- rnorm(1, 0, sqrt(k2 / (tau_TRUE * P)))
#   mu_j_TRUE[[i]] <- rnorm(J, mu_TRUE[i], sqrt(k1_TRUE / (tau_TRUE * P)))
#   y_final <- rnorm(n, mu_j_TRUE[[i]][x], sqrt(1 / tau_TRUE))
#   y<- y_final + y
# }

# now trying without splitting it
mu_TRUE <- rnorm(1, 0, sqrt(k2 / (tau_TRUE)))
mu_j_TRUE <- rnorm(J, mu_TRUE, sqrt(k1_TRUE / (tau_TRUE)))
y <- rnorm(n, mu_j_TRUE[x], sqrt(1 / tau_TRUE))

# Plot the data
p1 <- qplot(x, y) + geom_boxplot()
p1


# An alternative way to simulate the data - the code below should work as well (but not produce identical results) whether you use y here or the y above

# Fitting algorithm --------------------------------------------------------
n_iter <- 200
# Storing objects
# One mu_j for each group

# Some helpful values for marginal simulation later
ones_n <- rep(1, n)
M <- model.matrix(~ x - 1)
MtM <- tcrossprod(x = M) # Note tcrossprod(X) computes X%*%t(X)
tMM <- crossprod(x = M) # crossprod(X) computes t(X)%*%X

res <- rep(0, n)
mu_j <- matrix(0, nrow = J, ncol = P)
mu <- rep(0, P)
mu_j_store <- array(NA, dim = c(n_iter, J, P))
mu_store <- matrix(NA, ncol = P, nrow = n_iter)
k1_store <- rep(NA, n_iter)
tau_store <- rep(NA, n_iter)

# Starting values ---------------------------------------------------------
# Create values based on current value of k1
k1 <- 3
# k2 is given as the true one
Psi <- k1 * MtM + diag(n)
W_1 <- Psi + k2 * tcrossprod(x = ones_n)
p2_mu <- (t(ones_n) %*% solve(Psi, ones_n) + (P / k2))

df_y <- data.frame(y = y, group = as.numeric(x))
tau <- rgamma(1, 1 / alpha, beta)

# curve(dnorm(x,
#             mean = (p1_mu / p2_mu),
#             sd = sqrt((1/(tau * p2_mu)))), from = -2, to =  2)

# Gibbs sampling ----------------------------------------------------------
for (i in 1:n_iter) {
  if (i %% 10 == 0) print(i)

  # Sampling node/tree (stumps) parameters
  for (l in 1:P) {

    # Create the current set of fitted values and partial residuals
    res <- y - rowSums(mu_j[x, -l])

    # Sample a mu
    p1_mu <- t(ones_n) %*% solve(Psi, res)
    p2_mu <- (t(ones_n) %*% solve(Psi, ones_n) + (P / k2))
    mu[l] <- rnorm(1,
      mean = (p1_mu / p2_mu),
      sd = sqrt((1 / (tau * p2_mu)))
    )

    # Sample the mu_js for each group in each tree
    for (j in 1:J) {
      p1_muj <- (P * mu[l] / k1) + sum(res[x == j])
      p2_muj <- (n_j[j] + P / k1)
      mu_j[j, l] <- rnorm(1,
        mean = p1_muj / p2_muj,
        sd = sqrt((1 / (tau * p2_muj)))
      )
    }
  }

  # Update tau
  fits <- rowSums(mu_j[x, ])
  res_sq <- sum((y - fits)^2)
  mus_sq <- sum((mu_j - matrix(mu, nrow = J, ncol = P, byrow = TRUE))^2)
  mu_sq <- sum(mu^2)
  
  tau <- rgamma(
    1, (N + J*P + P) / 2 + alpha,
    0.5 *(res_sq + (P/k1)*mus_sq + (P/k2)*mu_sq) + beta
  )
  # tau <- tau_TRUE
  
  # Propose update for k1
  k1_new <- runif(1, 0, 20)
  Psi_new <- k1_new * MtM + diag(n)
  W_1_new <- Psi_new + k2 * tcrossprod(x = ones_n)
  p_new <- dmvnorm(y, sigma = W_1_new/tau, log = TRUE) + dunif(k1_new, 0, 20, log = TRUE)
  p_old <- dmvnorm(y, sigma = W_1/tau, log = TRUE) + dunif(k1, 0, 20, log = TRUE)
  a = exp(p_new - p_old)

  if(a > runif(1)) {
    k1 <- k1_new
    Psi <- Psi_new
    W_1 <- W_1_new
  }

  # Put all this into storage
  mu_store[i, ] <- mu
  mu_j_store[i, , ] <- mu_j
  k1_store[i] <- k1
  tau_store[i] <- tau
}

#stop()

# Check tau
qplot(tau_store) + geom_vline(xintercept = tau_TRUE)

# Check mu_j
qplot(apply(mu_j_store[n_iter,,], 1, 'sum'), mu_j_TRUE) + geom_abline()

# Check mu
qplot(rowSums(mu_store)) + geom_vline(xintercept = mu_TRUE)

# Check k1
qplot(k1_store) + geom_vline(xintercept = k1_TRUE)

# Check the fits
qplot(fits, y) + geom_abline()

# JAGS version ------------------------------------------------------------

# Model description -------------------------------------------------------
# We have y_{ij} ~ N(sum_k mu_jk, tau^-1)
# where i is obs and j is group, j=1,...,g, i=1,..,n_g
# k is the tree number k=1,...,m
# Let R_ijk = y_ij - sum_{l \ne k} mu_{jl}
# Then R_{ijk} ~ N(mu_jk, tau^-1)
# mu_jk ~ N(mu_k, k1 tau^{-1}/m)
# mu ~ N(0, k2 tau^{-1}/m)
# tau ~ Ga(alpha, beta)
# Hyper-parameter values --------------------------------------------------


# JAGS version ------------------------------------------------------------
model_code <- "
model
{
  # Likelihood
  for (i in 1:N) { # N is total number of obs
    y[i] ~ dnorm(fits[i], tau)
    fits[i] <- sum(mu_jk[group[i],1:P])
  }

  # Priors on means
  for (p in 1:P) { # P is the number of trees
    # Variance is k1/ (P * tau)
    # So precision is P * tau / k1
    for (j in 1:J) { # J is number of groups
      mu_jk[j, p] ~ dnorm(mu[p], P * tau / k1)
    }
    mu[p] ~ dnorm(0, P * tau / k2)
  }

  # Priors
  k1 ~ dunif(0, 20)
  tau ~ dgamma(alpha, beta)
}
"

# Set up the data
model_data <- list(
  N = N, P = 5,
  J = J, k2 = k2,
  group = x,
  y = y,
  alpha = alpha,
  beta = beta
)

# Choose the parameters to watch
model_parameters <- c("mu_jk", "mu", "k1", "tau", "fits")

# Run the model - can be slow
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code),
  n.chains = 1
)

# Check the fit
plot(model_run)

# Check tau
qplot(model_run$BUGSoutput$sims.list$tau) + 
  geom_vline(xintercept = tau_TRUE)

# Check mu_j
mu_j <- model_run$BUGSoutput$mean$mu_jk
qplot(rowSums(mu_j), mu_j_TRUE) + geom_abline()

# Check mu
mu <- model_run$BUGSoutput$sims.list$mu
qplot(rowSums(mu_store)) + geom_vline(xintercept = mu_TRUE)

# Check k1
k1 <- model_run$BUGSoutput$sims.list$k1
qplot(k1) + geom_vline(xintercept = k1_TRUE)

# Check the fits
fits <- model_run$BUGSoutput$mean$fits
qplot(fits, y) + geom_abline()
