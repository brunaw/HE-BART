library(tidyverse)
library(MASS)

det2 <- function(k_1_d, k_2_d, M_d) {
  n <- nrow(M_d)
  n_j <- colSums(M_d)
  tMM <- crossprod(x = M_d)
  Psi_tilde_inv <- diag(n) - M_d%*%diag(k_1_d/(1 + k_1_d*n_j))%*%t(M_d)
  return(log(k_2_d) + log(1/k_2_d + sum(Psi_tilde_inv)) + sum(log(1 + k_1_d*n_j)))
}

inv2 <- function(k_1, k_2, M) {
  n <- nrow(M)
  n_j <- colSums(M)
  Psi_tilde_inv <- diag(n) - M%*%diag(k_1/(1 + k_1*n_j))%*%t(M)
  k_3 <- 1/k_2 + sum(Psi_tilde_inv)
  Psi_row_sums <- rowSums(Psi_tilde_inv)
  Psi_col_sums <- colSums(Psi_tilde_inv)
  inv_M <- Psi_tilde_inv - tcrossprod(Psi_row_sums, Psi_col_sums) / k_3
  return(inv_M)
}

sim_y <- function(M, mu_mu = 0, taud = tau, k1d = k1, k2d = k2){
  n <- nrow(M)
  ones_n <- rep(1, n)
  MtM <- tcrossprod(x = M) # Note tcrossprod(X) computes X%*%t(X)
  tMM <- crossprod(x = M) # crossprod(X) computes t(X)%*%X
  W_1 <- k1d * MtM + k2d * tcrossprod(x = ones_n) + diag(n) 
  y <- mvrnorm(1, rep(mu_mu, n), W_1/taud)
  y
}

# Andrew's simulation function (depends on sim_y())
sim_a <- function(m, n = 1000, k1 = 8){
  alpha = 0.5; beta = 1; mu_mu = 0;  k2 = 10
  tau <- rgamma(1, 1/alpha, beta)
  alloc <- factor(sample(1:m, size = n, replace = TRUE))
  X1 <- runif(n)
  group_left  <- factor(alloc[X1 < 0.5])
  group_right <- factor(alloc[!(X1 < 0.5)])
  M_1_left    <- model.matrix(~ group_left - 1)
  M_1_right   <- model.matrix(~ group_right - 1)
  
  y_left <- sim_y(M_1_left, taud = tau, k1d = k1, k2d = k2)
  y_right <- sim_y(M_1_right, taud = tau, k1d = k1, k2d = k2)
  
  data <- data.frame(y = y_left, group = group_left, X1 = X1[X1 < 0.5]) %>% 
    bind_rows(data.frame(y = y_right, group = group_right, X1 = X1[!(X1 < 0.5)]))
  return(list(data = data, tau = tau)) 
}

# Bruna's simulation function
sim_b <- function(m, n = 1000, k1 = 8 ){
  alpha = 0.5; beta = 1; mu_mu = 0; k2 = 10
  alloc <- sample(1:m, size = n, replace = TRUE)
  X1 <- runif(n)
  tau <- rgamma(1, 1/alpha, beta) # 1.29
  mu <- rnorm(2, mu_mu, sqrt(k2/tau))
  
  muj_1 <- rnorm(m, mu[1], sqrt(k1/tau))
  muj_2 <- rnorm(m, mu[2], sqrt(k1/tau))
  
  y <- vector(length = n)
  for(i in 1:n) {
    curr_mean <- if(X1[i] < 0.5) { muj_1
    } else { muj_2 }
    y[i] <- rnorm(1, curr_mean[alloc[i]], sd = sqrt(1/tau))
  }
  data <- data.frame(X1, y, group = alloc) 
  return(list(data = data, tau = tau)) 
}

# Marginal distribution evalution function (with a stump)
ev_marg <- function(k_1, data, k_2 = 10, 
                    mu_mu = 0, alpha = 0.5, 
                    beta = 1){
  M <- model.matrix(~ factor(data$group) - 1)
  y <-  data$y
  n <- nrow(M)
  W_0 <- rep(mu_mu, n)
  ymW_0 <- y - W_0
  
  term_1 <- -(n/2)*log(2*pi)
  term_2 <- - 0.5 * det2(k_1_d = k_1, k_2_d = k_2, M_d = M)
  term_3 <- lgamma(n/2 + alpha)
  
  term_4 <- - (n/2 + alpha)*log(0.5 * t(ymW_0)%*%inv2(k_1, k_2, M)%*%ymW_0 + beta)
  all <- term_1 + term_3 + term_2 + term_4
  data.frame(term_1, term_2, term_3, term_4, all, k_1)
}

# True k1 = 8 
set.seed(2021)
data_a <- sim_a(m = 10, n = 1000)
data_b <- sim_b(m = 10, n = 1000)

k1_seq <- seq(0.01, 25, by = 0.1)
data_a <- data_a$data %>% 
  mutate(node = ifelse(X1 < 0.5, "l", "r")) %>% 
  split(.$node)

data_b <- data_b$data %>% 
  mutate(node = ifelse(X1 < 0.5, "l", "r")) %>% 
  split(.$node)
# Running the functional for a set of k1s and plotting 
ev_margs_a_l <- map_dfr(k1_seq, ev_marg, data = data_a$l)
ev_margs_a_r <- map_dfr(k1_seq, ev_marg, data = data_a$r)
ev_margs_b_l <- map_dfr(k1_seq, ev_marg, data = data_b$l)
ev_margs_b_r <- map_dfr(k1_seq, ev_marg, data = data_b$r)

ev_margs_a_l %>% mutate(node = "l", type = "Andrew") %>% 
  bind_rows(ev_margs_a_r %>% mutate(node = "r", type = "Andrew")) %>% 
  bind_rows(ev_margs_b_l %>% mutate(node = "l", type = "Bruna") %>% 
              bind_rows(ev_margs_b_r %>% mutate(node = "r", type = "Bruna"))) %>% 
  gather(key, value, -k_1, -type, -node) %>%
  group_by(k_1, type, key) %>% 
  mutate(value = sum(value)) %>% 
  
  filter(!str_detect(key, "1|3")) %>% 
  group_by(type, key) %>% 
  mutate(m = max(value), 
         k_1_max = ifelse(value == m, k_1, NA)) %>% 
  ggplot(aes(k_1, value)) +
  facet_wrap(type~key, scales = 'free', ncol = 3) +
  geom_point() +
  geom_hline(aes(yintercept = m), colour = "red") +
  geom_vline(aes(xintercept = k_1_max), colour = "red") +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  labs(title = "Term 2 = - 0.5 * det2(k_1_d = k_1, k_2_d = k_2, M_d = M) ; 
Term 4 = - (n/2 + alpha)*log(0.5 * t(ymW_0)%*%inv2(k_1, k_2, M)%*%ymW_0 + beta);
All = sum of all terms") +
  theme_bw(12) 

ggsave(file = "results/ll_k1_terms.png", 
       height = 6, width = 10)


