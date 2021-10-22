library(MASS)
library(tidyverse)
set.seed(2022)
alpha = 0.5; beta = 1; mu_mu = 0; k1 = 8 ; k2 = 10
tau <- rgamma(1, 1/alpha, beta) #  1.043341

n <- 1e3 # You can increase this, but some of the comparisons below get very slow
m <- 10
alloc <- factor(sample(1:m, size = n, replace = TRUE))
X1 <- runif(n)

group_left  <- factor(alloc[X1 < 0.5])
group_right <- factor(alloc[!(X1 < 0.5)])
M_1_left    <- model.matrix(~ group_left - 1)
M_1_right   <- model.matrix(~ group_right - 1)

sim_y <- function(M, mu_mu = 0, taud = tau, 
                  k1d = k1, k2d = k2){
  n <- nrow(M)
  ones_n <- rep(1, n)
  MtM <- tcrossprod(x = M) # Note tcrossprod(X) computes X%*%t(X)
  tMM <- crossprod(x = M) # crossprod(X) computes t(X)%*%X
  W_1 <- k1d * MtM + k2d * tcrossprod(x = ones_n) + diag(n) 
  y <- mvrnorm(1, rep(mu_mu, n), W_1/taud)
  y
}

y_left <- sim_y(M_1_left)
y_right <- sim_y(M_1_right)

data <- data.frame(y = y_left, group = group_left, 
           X1 = X1[X1 < 0.5]) %>% 
  bind_rows(
    data.frame(y = y_right, group = group_right, 
               X1 = X1[!(X1 < 0.5)]))

data %>% 
  group_by(group) %>% 
  summarise(m = mean(y))

iter <- 1000
lim <- 200
formula <- y ~ X1
group_variable <- "group"
pars <- list(
  k1 = 6, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)

# Andrew's simulated data 
# Correctly estimates k_1
# 2-node fixed tree ------- 
m0_k1_samp_0 <- bcart_fixed_k1(
  formula, 
  dataset = data,
  iter = iter, 
  group_variable, pars, 
  scale = FALSE, 
  min_u = 0, max_u = 15)

current_tree = m0_k1_samp_0$final_tree
m0_k1_samp_0$tau_post %>% mean()
m0_k1_samp_0$sampled_k1 %>% summary()
plot(density(m0_k1_samp_0$sampled_k1))

res_sample_k1 <- data.frame(
  taus = m0_k1_samp_0$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m0_k1_samp_0$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m0_k1_samp_0$mus, ~{.x[2]})[-c(1:lim)],
  k1 = m0_k1_samp_0$sampled_k1[-c(2:lim)]
)

res_sample_k1 %>% 
  group_by(k1) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(mean_k1 = mean(k1)) %>% 
  ggplot(aes(x = k1)) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean_k1), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 8), 
             colour = "blue", size = 1) +
  theme_bw()

ggsave(file = "results/k1_samples_marg.png")

# res_sample_k1 %>% 
#   # group_by(k1) %>%
#   # slice(1) %>%
#   ungroup() %>% 
#   mutate(mean_k1 = mean(k1), iteration = 1:n()) %>% 
#   ggplot(aes(y = k1, x = iteration)) +
#   geom_point(fill = "grey90") +
#   geom_hline(aes(yintercept = mean_k1), 
#              colour = "red", linetype = "dashed", size = 1) +
#   theme_bw()
# ggsave(file = "results/k1_samples_chain_big.png")


res_sample_k1 %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/tau_samp_k1_marg.png")


k1_tau <- 8/tau
res_sample_k1 %>% 
  group_by(k1) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(mean_k1_tau = mean(k1/taus)) %>% 
  ggplot(aes(x = sqrt(k1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = sqrt(mean_k1_tau)), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(k1_tau)), 
             colour = "blue", size = 1) +
  #xlim(0, 2) +
  theme_bw()
ggsave(file = "results/k1_tau_marg.png")


# mu = c(-1.628563, 1.483883)

res_sample_k1 %>% 
  select(2:3) %>% 
  gather() %>% 
  #mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key) + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  # geom_vline(aes(xintercept = true_mean), 
  #            colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/mu_samp_k1_marg.png")
#--------------------------------------------------------------------
# # My simulated data 
# # Correctly estimates k_1? 
iter = 1000
# 2-node fixed tree -------
m0_k1_bruna <- bcart_fixed_k1(
  formula,
  dataset = data_k1,
  iter = iter,
  group_variable, pars,
  scale = FALSE,
  min_u = 0, max_u = 15)

n_distinct(data_k1$group)
mean(m0_k1_bruna$sampled_k1)
 
