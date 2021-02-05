psi <- (k + diag(N))
det_psi <- det(psi)
in_psi <- solve(psi)

params <- list(
  N = N, mu_j  = mu_j, y = y, 
  group = group, beta = beta, 
  alpha = alpha, tau = tau, 
  tau_mu = tau_mu, 
  mu_mu = mu_mu,
  m = M, 
  nj = nj,
  mu_main = mu_main, k = k,
  det_psi = det_psi,
  in_psi = in_psi,
  psi = psi
)

cond_y <- function(params){
  N <-  length(y)
  k <- params$k
  det_psi <- params$det_psi
  psi <- params$psi
  in_psi <- params$in_psi
  beta <- params$beta
  alpha <- params$alpha
  #group <- params$group
  #mu_j <- params$mu_j
  #tau <- params$tau
  #tau_mu <- params$tau_mu
  mu_mu <- params$mu_mu
  y <- params$y
  m <-  params$m
  
  if(k != 0.5){
    psi <- (k + diag(N))
    det_psi <- det(psi)
    in_psi <- solve(psi)
  }

  M_mat <- model.matrix(y ~ factor(group) - 1)
  
  mm <- rep(1, N) %*% psi %*% rep(1, N)
  multiplier <- (det_psi^(-1/2)) * (k^(1/2)) * (1/(k + mm^(1/2)))
  
  inner <-  ((k* mu_mu^2) + (y %*% in_psi %*% y) + 
               (k*mu_mu + y %*% in_psi %*% M_mat %*% rep(1, m))^2/(k + mm))
  
  # Missing multiplier? that's likely 
  # turn to log 
  
  p_cond_y <- log(multiplier) + (log( beta * (1/2)*(inner)))*(-(N/2 + alpha))
  
  #p_exp <- exp(p_cond_y)
  return(p_cond_y)
}

cond_y(params)

# Profiling parameters -----------------------------
df_dens <- data.frame(
  k = seq(0, 10, length.out = 25), 
  mu = seq(-3, 3, length.out = 25), 
  beta = seq(0, 5, length.out = 25)) %>% 
  as_tibble() %>% 
  dplyr::mutate(
    params_v = list(params), 
    params_v = map2(params_v, mu, 
                    ~list_modify(.x, mu_mu = .y)),
    params_v2 = map2(params_v, k, 
                    ~list_modify(.x, k = .y)), 
    params_v3 = map2(params_v, beta, 
                     ~list_modify(.x, beta = .y))) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(cond_y_mu = cond_y(params_v), 
         cond_y_k = cond_y(params_v2),
         cond_y_beta = cond_y(params_v3))


df_dens %>% 
  ggplot(aes(x = mu, y = unlist(cond_y_mu))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = mu_mu, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = expression(mu[mu]), y = "log-likelihood", 
       title = expression("Conditional distribution for "~y))
ggsave(file = "book2/img/ll_prof_mu.png")

df_dens %>% 
  ggplot(aes(x = k, y = unlist(cond_y_k))) +
  geom_line() +
  theme_light(18) +
  geom_vline(xintercept = k, linetype = 2, 
             colour = 'tomato', size = 0.65) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = expression(k), y = "log-likelihood", 
       title = expression("Conditional distribution for "~y))
ggsave(file = "book2/img/ll_prof_k.png")

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
ggsave(file = "book2/img/ll_prof_beta.png")


