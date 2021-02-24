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

nj = sample(50:300, M, replace = TRUE) # Set the number of obs in each group between 5 and 10
N = sum(nj)
group = rep(1:M, times = nj)
mu = rep(mu_j, nj)
#M1 = model.matrix(mu ~ factor(group) - 1)
y = rnorm(N, mean = mu, sd = sqrt(1/tau))

M_mat <- model.matrix(mu ~ factor(group) - 1)
psi <- (k1 * M_mat %*% t(M_mat)) + diag(N)
det_psi <- det(psi)
in_psi <- solve(psi)
W_1 <- (k2 * diag(x = 1, nrow = N, ncol = 1) %*% t(diag(x =1, nrow = N, ncol = 1))) + psi
in_W_1 <- solve(W_1)


params <- list(
  N = N, mu_j  = mu_j, y = y, 
  group = group, beta = beta, 
  alpha = alpha, tau = tau, 
  tau_mu = tau_mu, 
  mu_mu = mu_mu,
  m = M, 
  nj = nj,
  mu_main = mu_main, 
  k1 = k1,
  k2 = k2, 
  W_1 = W_1, 
  in_W_1 = in_W_1, 
  det_psi = det_psi,
  in_psi = in_psi,
  psi = psi
)

cond_y <- function(params){
  N <-  length(y)
  k1 <- params$k1
  k2 <- params$k2
  det_psi <- params$det_psi
  psi <- params$psi
  in_psi <- params$in_psi
  beta <- params$beta
  alpha <- params$alpha
  mu_mu <- params$mu_mu
  mu_main <- params$mu_main
  y <- params$y
  m <-  params$m
  W_1 <-  params$W_1
  in_W_1 <-  params$in_W_1
  
  M_mat <- model.matrix(y ~ factor(group) - 1)
  
  if(k1 != 0.5){
    psi <- (k1 * M_mat %*% t(M_mat)) + diag(N)
    det_psi <- det(psi)
    in_psi <- solve(psi)
    W_1 <- (k2 * diag(x = 1, nrow = N, ncol = 1) %*% t(diag(x =1, nrow = N, ncol = 1))) + psi
    in_W_1 <- solve(W_1)
  }
  
  W_0 <- rep(mu_main, N)

  inner <- (t((y - W_0)) %*% in_W_1 %*% (y - W_0)) + beta
  #  all <-  (det(W_1)^(-1/2) * inner^(-(N/2 + alpha)))

  p_cond_y <- (-1/2) * log(det(W_1)) + (log(inner)*(-(N/2 + alpha))) +
    lgamma(N/2 + alpha)
  
  #p_exp <- exp(p_cond_y)
  return(p_cond_y)
}

# params$alpha <- 0.8
# cond_y(params)

# Profiling parameters -----------------------------
df_dens <- data.frame(
  k1 = seq(0, 10, length.out = 25), 
  mu = seq(-3, 3, length.out = 25), 
  beta = seq(0, 5, length.out = 25),
  alpha = seq(0, 5, length.out = 25)) %>% 
  as_tibble() %>% 
  dplyr::mutate(
    params_v = list(params), 
    params_v = map2(params_v, mu, 
                    ~list_modify(.x, mu_main = .y)),
    params_v2 = map2(params_v, k1, 
                    ~list_modify(.x, k1 = .y)), 
    params_v3 = map2(params_v, beta, 
                     ~list_modify(.x, beta = .y)),
    params_v4 = map2(params_v, alpha, 
                     ~list_modify(.x, alpha = .y))) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(cond_y_mu = cond_y(params_v), 
         cond_y_k1 = cond_y(params_v2),
         cond_y_beta = cond_y(params_v3),
         cond_y_alpha = cond_y(params_v4))


df_dens %>% 
  ggplot(aes(x = mu, y = unlist(cond_y_mu))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = mu_mu, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = expression(mu), y = "log-likelihood", 
       title = expression("Conditional distribution for "~y))
#ggsave(file = "book2/img/ll_prof_mu.png")

df_dens %>% 
  ggplot(aes(x = k1, y = unlist(cond_y_k1))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = k1, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = expression(k), y = "log-likelihood", 
       title = expression("Conditional distribution for "~y))
#ggsave(file = "book2/img/ll_prof_k.png")

df_dens %>% 
  ggplot(aes(x = beta, y = unlist(cond_y_beta))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = beta, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = expression(beta), y = "log-likelihood", 
       title = expression("Conditional distribution for "~y))
#ggsave(file = "book2/img/ll_prof_beta.png")


df_dens %>% 
  ggplot(aes(x = alpha, y = unlist(cond_y_alpha))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = alpha, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = expression(alpha), y = "log-likelihood", 
       title = expression("Conditional distribution for "~y))
#ggsave(file = "book2/img/ll_prof_alpha.png")


