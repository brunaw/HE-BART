library(tidyverse)
library(reshape2)
library(mvtnorm)
library(gridExtra)
set.seed(123)

# Model description -------------------------------------------------------
# We have R_{ij} ~ N(mu_j, tau^{-1})
# mu_j ~ N(mu, k1 tau^{-1}/P)
# mu ~ N(0, k2 tau^{-1}/{})
# tau ~ Ga(alpha, beta)
# Hyper-parameter values --------------------------------------------------
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
M <- model.matrix(~ x - 1)
# Mean shrinkage parameter (always fixed)
k2 <- 5
# Conditional simulation
k1_TRUE <- 8
tau_TRUE <- rgamma(1, alpha, beta)
y<- vector(length = n)
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
mu_TRUE <- rnorm(1, 0, sqrt(k2 / (tau_TRUE )))
mu_j_TRUE <- rnorm(J, mu_TRUE, sqrt(k1_TRUE / (tau_TRUE)))
y <- rnorm(n, mu_j_TRUE[x], sqrt(1 / tau_TRUE))

# Plot the data
p1 <- qplot(x, y) + geom_boxplot()
p1

# Some helpful values for marginal simulation later
ones_n <- rep(1, n)
MtM <- tcrossprod(x = M) # Note tcrossprod(X) computes X%*%t(X)
tMM <- crossprod(x = M) # crossprod(X) computes t(X)%*%X

# An alternative way to simulate the data - the code below should work as well (but not produce identical results) whether you use y here or the y above

# Fitting algorithm --------------------------------------------------------
n_iter <- 50
# Storing objects
# One mu_j for each group
mu_j <- rep(NA, J)
mu_store <- data.frame(
  tree = rep(1:P, n_iter),
  iter = rep(1:n_iter, each = 5)
) %>% 
  as_tibble() %>% 
  mutate(mu_sample = NA)

mu_j_store <- data.frame(
  tree = rep(1:P, n_iter),
  iter = rep(1:n_iter, each = 5)
) %>% 
  as_tibble() %>% 
  mutate(mu_j_sample = list(vector(length = J)))

k1_store <- rep(NA, n_iter)
tau_store <- rep(NA, n_iter)

# Starting values ---------------------------------------------------------
# Create values based on current value of k1
k1 <- 3
# k2 is given as the true one
Psi <- k1 * MtM + diag(n)
W_1 <- Psi + k2 * tcrossprod(x = ones_n) 
p2_mu <- (t(ones_n) %*% solve(Psi, ones_n) + (P/ k2))

df_y <- data.frame(y = y, group = as.numeric(x))
tau <- rgamma(1, 1/alpha, beta)

# curve(dnorm(x,
#             mean = (p1_mu / p2_mu),
#             sd = sqrt((1/(tau * p2_mu)))), from = -2, to =  2)

# Gibbs sampling ----------------------------------------------------------
for (i in 1:n_iter) {
  if (i %% 10 == 0) print(i)
  
  # Defines residuals to be used ------------
  df_res <- data.frame(res = y, group = as.numeric(x))
  
  # Sampling node/tree (stumps) parameters
  #---------------------------------------------
  for(l in 1:P){
    res <- df_res$res
    p1_mu <- t(ones_n) %*% solve(Psi, res) 
    
    # Sample a mu for each tree
    # check that these parameters match to the maths
    mu <- rnorm(1,
                mean = (p1_mu / p2_mu),
                sd = sqrt((1/(tau * p2_mu)))
    )
    #mu <- mu/P
    #mu <- mu_TRUE/P
    
    # Stores it 
    mu_store <- mu_store %>% 
      mutate(mu_sample = 
               ifelse(tree == l & iter == i, mu, 
                      mu_sample))
    
    # Sample the mu_js for each tree
    for (j in 1:J) {
      p1_muj <- (P*mu / k1) + (mean(res[x == j]) * n_j[j])
      p2_muj <- (n_j[j] + P/k1)
      mu_j[j] <- rnorm(1,
                       mean = p1_muj / p2_muj,
                       sd = sqrt((1/(tau* p2_muj)))
      )
    }
    #mu_j <- mu_j/P
    # 
    #mu_j <- mu_j_TRUE/P
    
    # Stores it
    mu_j_store <- mu_j_store %>% 
      mutate(mu_j_sample = 
               ifelse(tree == l & iter == i, list(mu_j), mu_j_sample))
    
    # calculating partial residuals for next tree
    df_res$res <- mu_j_store %>% 
      filter(iter == i, tree == l) %>% 
      select(-iter) %>% 
      unnest(mu_j_sample) %>% 
      mutate(group = 1:J) %>% 
      left_join(df_res, by = "group") %>% 
      mutate(new_res = res - mu_j_sample) %>% 
      pull(new_res)
  }
  #---------------------------------------------
  
  f_is <- mu_j_store %>%
    filter(iter == i) %>%
    select(-iter) %>%
    unnest(mu_j_sample) %>%
    mutate(group = rep(1:J, P)) %>%
    group_by(group) %>%
    summarise(sum_mu = sum(mu_j_sample))

  
  # Sample tau
  res_sq <- f_is %>% 
    left_join(df_y, by = "group") %>% 
    summarise(sum_sq = sum((y - sum_mu)^2)) %>% 
    pull(sum_sq) 
  
  # sum sq of mu and mu_js 
  mus_sq <- mu_j_store %>% 
    filter(iter == i) %>% 
    select(-iter) %>% 
    unnest(mu_j_sample) %>% 
    left_join(
      mu_store %>% 
        filter(iter == i) %>% 
        select(-iter), by = "tree"
    ) %>% 
    summarise(sum_sq = sum((mu_j_sample - mu_sample)^2)) %>% 
    pull(sum_sq)
  
  mu_sq <- mu_store %>% 
    filter(iter == i) %>% 
    select(-iter) %>% 
    summarise(sum_sq = sum(mu_sample^2)) %>% 
    pull(sum_sq)
  
  #---------------------------------------------
  # Tau and k1 are only updated after we have the P trees 
  tau <- rgamma(
    1, (N + J*P + P) / 2 + alpha,
    0.5 *(res_sq + (P/k1)*mus_sq + (P/k2)*mu_sq) + beta
  )
  
  
  # Propose update for k1
  # k1_new <- runif(1, 0, 20)
  # Psi_new <- k1_new * MtM + diag(n)
  # W_1_new <- Psi_new + k2 * tcrossprod(x = ones_n)
  # p_new <- dmvnorm(y, sigma = W_1_new/tau, log = TRUE) + dunif(k1_new, 0, 20, log = TRUE)
  # p_old <- dmvnorm(y, sigma = W_1/tau, log = TRUE) + dunif(k1, 0, 20, log = TRUE)
  # a = exp(p_new - p_old)
  # 
  # if(a > runif(1)) {
  #   k1 <- k1_new
  #   Psi <- Psi_new
  #   W_1 <- W_1_new
  # }
  
  
  # Put all this into storage
  k1_store[i] <- k1
  tau_store[i] <- tau
}


#----------------------------------------------------------------------------
# Case 1: ok, tau is sampled correctly
# Case 2: ok, tau is sampled correctly but with a heavy tail; tau underestimated
# Case 3 (mu sampled, mu_j fixed): tau is mostly underestimated/overestimated
# Case 3 with mu/p: tau correctly sampled
# Case 2 with mu_j/p: tau very underestimated
# mu_js are slightly off

# There is a problem with both mu and mu_j

# mu_j is being sampled correctly
# Plot together
df_k1 <- data.frame(R = k1_store) %>% 
  group_by(R) %>% 
  slice(1)

p1 <- ggplot(df_k1, aes(
  x = R
)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = k1_TRUE) +
  ggtitle("k1 posterior")
p1

df_tau <- data.frame(R = tau_store)

p2 <- ggplot(df_tau, aes(
  x = R
)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = tau_TRUE) +
  ggtitle("tau posterior")
p2
#--------------------------------------------------------
mu_true_df <- data.frame(mu_TRUE, tree = 1)

mu_store_2 <- mu_store %>% 
  group_by(iter) %>%
  summarise(m = sum(mu_sample)) 
  #left_join(mu_true_df, by = 'tree')

p3 <- ggplot(mu_store_2, aes(
  x = m
  )) +
  geom_density(alpha = 0.4) +
  geom_vline(data = mu_true_df, aes(xintercept = mu_TRUE)) +
  ggtitle("Mu posterior")
p3
#--------------------------------------------------------------------
# # Finally mu_j
df_mu_j <- mu_j_store %>% 
  unnest(mu_j_sample) %>% 
  mutate(group = rep(1:J,  n_iter * P)) %>%
  group_by(iter, group) %>% 
  summarise(value = sum(mu_j_sample))

mu_j_truth <- data.frame(
  group = 1:J,
  value = mu_j_TRUE
)

p4 <- ggplot(df_mu_j, aes(x = value, fill = factor(group))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~group, scales = "free") +
  ggtitle("Mu_j posteriors") +
  geom_vline(
    data = mu_j_truth,
    aes(xintercept = value),
    col = "red"
  ) + 
  guides(fill = FALSE)
p4
#---------------------------------------------------------------
library(patchwork)

p1 + p2 + plot_layout(nrow = 2)
p3 + p4 + plot_layout(nrow = 2)

# JAGS code is not complete 
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

library(R2jags)
# JAGS version ------------------------------------------------------------
model_code <- "
model
{
  # Likelihood
  for (i in 1:N) { # N is total number of obs
    y[i] ~ dnorm(sum(mu_jk[group[i],1:P]), tau)
  }
  
  # Priors on means
  for (p in 1:P) { # P is the number of trees
    # Variance is k1/ (G * tau)
    # So precision is G * tau / k1
    for (j in 1:J) { # J is number of groups
      mu_jk[j, p] ~ dnorm(mu[p], G * tau / k1)
    }
    mu[p] ~ dnorm(0, G * tau / k2)
    # Compute the partial residuals too
    # R[1:N, p] <- y - 
  }

  # Priors
  k1 ~ dunif(0, 20)
  tau ~ dgamma(alpha, beta)
}
"

# Set up the data
model_data <- list(N = N, P = 5, 
                   J = J, k2 = k2,
                   group = x,
                   y = y,
                   alpha = alpha,
                   beta = beta)

# Choose the parameters to watch
model_parameters <- c( "mu_jk", "mu","k1", "tau")

# Run the model - can be slow
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code),
  n.chains = 1)
​
plot(model_run)
​
# Set up storage ----------------------------------------------------------
​
# Check it works - compare the mu values
mu_est <- model_run$BUGSoutput$mean$mu
sum(mu_est)
mu_TRUE
​
# Compare tau
model_run$BUGSoutput$mean$tau
tau_TRUE
​
# Compare mu_j
mu_j_est <- model_run$BUGSoutput$mean$mu_jk
apply(mu_j_est, 1, 'sum')
mu_j_TRUE
df_mu_j %>% 
  filter(iter == n_iter) %>% 
  left_join(df_y, by = 'group') %>% 
  summarise(ssq = sum((y - value)^2))

df_y %>% 
  group_by(group) %>% 
  mutate(m = mean(y)) %>% 
  ungroup() %>% 
  summarise(ssq = sum((y - m)^2))
# --------------------------------------------------------------------
