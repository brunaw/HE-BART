# Simulate some data from a stump with a set form, then run some Gibbs sampling on it to get the posterior estimates of the parameters
# Load in packages --------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(reshape2)
library(R2jags)
library(mvtnorm)
library(gridExtra)
set.seed(123)

# Model description -------------------------------------------------------
# We have y_{ij} ~ N(mu_j, tau^{-1})
# mu_j ~ N(mu, k1 tau^{-1})
# mu ~ N(0, k2 tau^{-1})
# tau ~ Ga(alpha, beta)
# Hyper-parameter values --------------------------------------------------
# Number of obs
n <- 1e3 # You can increase this, but some of the comparisons below get very slow
m <- 10 # Number of groups

# Other hyper-parameters defined above
alpha <- 1
beta <- 1

# Simulate data -----------------------------------------------------------
# Create a factor with the groups
x <- factor(sample(1:m, size = n, replace = TRUE))
# Group sizes
n_j <- table(x)
# Group allocation matrix
M <- model.matrix(~ x - 1)
# Mean shrinkage parameter (always fixed)
k2 <- 8
# Conditional simulation
k1_TRUE <- 8
tau_TRUE <- rgamma(1, alpha, beta)
mu_TRUE <- rnorm(1, 0, sqrt(k2 / tau_TRUE))
mu_j_TRUE <- rnorm(m, mu_TRUE, sqrt(k1_TRUE / tau_TRUE))
y <- rnorm(n, mu_j_TRUE[x], sqrt(1 / tau_TRUE))

# Plot the data
p1 <- qplot(x, y) + geom_boxplot()

# Some helpful values for marginal simulation later
ones_n <- rep(1, n)
MtM <- tcrossprod(x = M) # Note tcrossprod(X) computes X%*%t(X)
tMM <- crossprod(x = M) # crossprod(X) computes t(X)%*%X

# An alternative way to simulate the data - the code below should work as well (but not produce identical results) whether you use y here or the y above
# y <- mvrnorm(1, rep(mu_mu, n), W_1/tau)
# p2 <- qplot(x, y2) + geom_boxplot()

# Set up storage ----------------------------------------------------------

n_iter <- 500
mu_j <- rep(NA, m)
mu_store <- rep(NA, n_iter)
k1_store <- rep(NA, n_iter)
mu_j_store <- matrix(NA, ncol = m, nrow = n_iter)
tau_store <- rep(NA, n_iter)
# Starting values ---------------------------------------------------------
# Create values based on current value of k1
k1 <- 1
Psi <- k1 * MtM + diag(n)
W_1 <- Psi + k2 * tcrossprod(x = ones_n)
p1_mu <- t(ones_n) %*% solve(Psi, y)
p2_mu <- t(ones_n) %*% solve(Psi, ones_n) + 1 / k2

# Gibbs sampling ----------------------------------------------------------
for (i in 1:n_iter) {
  if (i %% 10 == 0) print(i)
  
  # Sample tau
  tau <- rgamma(
    1, n / 2 + alpha,
    0.5 * t(y) %*% solve(W_1, y) + beta
  )
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
    p1_mu <- t(ones_n) %*% solve(Psi, y)
    p2_mu <- t(ones_n) %*% solve(Psi, ones_n) + 1 / k2
  }
  
  # Sample mu
  mu <- rnorm(1,
              mean = p1_mu / p2_mu,
              sd = sqrt(1 / tau * 1 / p2_mu)
  )
  
  # Sample the mu_j
  for (j in 1:m) {
    p1_muj <- mu / k1_TRUE + mean(y[x == j]) * n_j[j]
    p2_muj <- n_j[j] + 1 / k1_TRUE
    mu_j[j] <- rnorm(1,
                     mean = p1_muj / p2_muj,
                     sd = sqrt(1 / (tau * p2_muj))
    )
  }
  
  # Put all this into storage
  k1_store[i] <- k1
  tau_store[i] <- tau
  mu_store[i] <- mu
  mu_j_store[i, ] <- mu_j
}

# JAGS version ------------------------------------------------------------

model_code <- "
model {
 for (i in 1:N) {
   y[i] ~ dnorm(mu_j[group[i]], tau)
 }
 for (j in 1:M) {
   mu_j[j] ~ dnorm(mu, tau/k1)
 }
 mu ~ dnorm(0, tau/k2)
 tau ~ dgamma(alpha, beta)
 k1 ~ dunif(0, 20)
}
"

model_data <- list(
  y = y,
  N = n,
  M = m,
  group = x,
  k2 = k2,
  alpha = alpha,
  beta = beta
)

model_run <- jags(
  data = model_data,
  parameters.to.save = c("mu", "mu_j", "tau", "k1"),
  model.file = textConnection(model_code)
)

# plot(model_run)
# Compare JAGS and R ------------------------------------------------------
k1_store_jags <- model_run$BUGSoutput$sims.list$k1
tau_store_jags <- model_run$BUGSoutput$sims.list$tau
mu_store_jags <- model_run$BUGSoutput$sims.list$mu
mu_j_store_jags <- model_run$BUGSoutput$sims.list$mu_j

# Plot together
df_k1 <- melt(data.frame(
  JAGS = model_run$BUGSoutput$sims.list$k1,
  R = k1_store
))
ggplot(df_k1, aes(
  x = value,
  fill = variable
)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = k1_TRUE) +
  ggtitle("k1 posteriors")

df_tau <- melt(data.frame(
  JAGS = model_run$BUGSoutput$sims.list$tau,
  R = tau_store
))
ggplot(df_tau, aes(
  x = value,
  fill = variable
)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = tau_TRUE) +
  ggtitle("Tau posteriors")

df_mu <- melt(data.frame(
  JAGS = model_run$BUGSoutput$sims.list$mu,
  R = mu_store
))
ggplot(df_mu, aes(
  x = value,
  fill = variable
)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = mu_TRUE) +
  ggtitle("Mu posteriors")

# Finally mu_j
mu_j_JAGS <- melt(model_run$BUGSoutput$sims.list$mu_j)
mu_j_JAGS$Type <- "JAGS"
mu_j_R <- melt(mu_j_store)
mu_j_R$Type <- "R"
df_mu_j <- rbind(mu_j_JAGS, mu_j_R)
mu_j_truth <- data.frame(
  Var2 = 1:m,
  value = mu_j_TRUE
)

ggplot(df_mu_j, aes(x = value, fill = Type)) +
  geom_density(alpha = 0.4) +
  geom_vline(
    data = mu_j_truth,
    aes(xintercept = value),
    col = "red"
  ) +
  facet_wrap(~Var2, scales = "free") +
  ggtitle("Mu_j posteriors")