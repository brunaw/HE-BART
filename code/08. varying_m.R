#------------------------------------------------------
# This code creates datasets with different m values
# (number of groups) and checks how the models estimake k_1
# August,  2021
#------------------------------------------------------
library(MASS)
library(tidyverse)
source("code/00. simulation_functions.R")
files <- list.files("code/hmb/R") %>% 
  paste0("code/hmb/R/", .)
map(files, source)

source("code/mixedbart/R/bcart_fixed_k1.R")

#------------------------------------------------------
# Auxiliar functions
simulated_and_fit <- function(m, type = "a", 
                              iter = 300, 
                              formula = y ~ X1,
                              group_variable = "group", 
                              stump_d = FALSE){
  if(type == "a"){ data <- sim_a(m = m)
  } else if(type == "b"){ 
    data <- sim_b(m = m)
  } else if(type == 'stump'){
    data <- sim_stump(m = m)
  }
  pars <- list(k1 = 8, k2 = 10, alpha = 0.5, beta = 1, mu_mu = 0)
  
  m0_stump <- bcart_fixed_k1(
    formula, 
    dataset = data$data,
    iter = iter, 
    group_variable, pars, 
    scale = FALSE, 
    min_u = 0, max_u = 15, prior_k1 = FALSE, stump = TRUE)
  
  m0 <- bcart_fixed_k1(
    formula, 
    dataset = data$data,
    iter = iter, 
    group_variable, pars, 
    scale = FALSE, 
    min_u = 0, max_u = 15, prior_k1 = FALSE, stump = FALSE)

  return(list(m0 = m0, m0_stump = m0_stump, tau = data$tau))
}

extract_results <- function(model){
  tau_est <- mean(model$m0$tau_post)
  tau <- model$tau
  k1_est <- mean(unique(model$m0$sampled_k1))
  return(list(tau_est = tau_est, tau = tau, k1_est = k1_est,
              tau_full = model$m0$tau_post,
              k1_full = model$m0$sampled_k1))
}

extract_results_stump <- function(model){
  tau_est <- mean(model$m0_stump$tau_post)
  tau <- model$tau
  k1_est <- mean(unique(model$m0_stump$sampled_k1))
  return(list(tau_est = tau_est, tau = tau, k1_est = k1_est,
              tau_full = model$m0_stump$tau_post,
              k1_full = model$m0_stump$sampled_k1))
}

m_groups <- c(3, 7, 5, 10, 15, 18, 25, 30, 50, 75, 100)
ss <- simulated_and_fit(5, "stump")

safe_f <- safely(simulated_and_fit)

# Now doing it with a stump and not a tree structure

run_all <- tibble(
  m = rep(m_groups, 3),
  type = rep(c("a", "b", "stump"), each = length(m_groups))) %>%
  mutate(model = map2(m, type, safe_f))

saveRDS(run_all, "results/varying_m.rds")

results <- run_all2 %>%
  #slice(-c(13, 26)) %>%
  mutate(results_tree = map(map(model, "result"), extract_results),
         results_stump = map(map(model, "result"), extract_results_stump))


# Uncertainty plots ------------------------------------------------------------
# results %>% 
#   filter(type == 'stump') %>% 
#   mutate(k1_est = map(results_tree, "k1_full")) %>% 
#   select(m, type, k1_est) %>% 
#   unnest(k1_est) %>% 
#   group_by_all() %>% 
#   slice(1) %>% 
#   ggplot(aes(factor(m), k1_est)) +
#   geom_boxplot(alpha = 0.5, fill = 'salmon') +
#   facet_wrap(~type) + 
#   geom_hline(yintercept = 8, colour = 'blue', linetype = 'dashed') +
#   #geom_point(aes(colour = type)) +
#   #scale_x_continuous(breaks = c(3, 7, 5, 10, 12, 15, 18, 20, 25, 30)) +
#   scale_y_continuous(breaks = 1:15) +
#   labs(y = "Estimated k1", x = "M (number of groups)") +
#   theme_bw()
# 
# ggsave(file = "results/varying_m.png", height = 7, width = 10)


results %>% 
  filter(type == 'stump') %>% 
  mutate(k1_est = map(results_stump, "k1_full")) %>% 
  select(m, type, k1_est) %>% 
  unnest(k1_est) %>% 
  group_by_all() %>% 
  slice(1) %>% 
  ggplot(aes(factor(m), k1_est)) +
  geom_boxplot(alpha = 0.5, fill = 'salmon') +
  facet_wrap(~type) + 
  geom_hline(yintercept = 8, colour = 'blue', linetype = 'dashed') +
  #geom_point(aes(colour = type)) +
  #scale_x_continuous(breaks = c(3, 7, 5, 10, 12, 15, 18, 20, 25, 30)) +
  scale_y_continuous(breaks = 1:15) +
  labs(y = "Estimated k1", x = "M (number of groups)") +
  theme_bw()

ggsave(file = "results/varying_m_stump.png", height = 7, width = 10)
# ------------------------------------------------------------------------------
# results %>% 
#   mutate(tau_est = map(results_tree, "tau_est"),
#          tau = map(results_tree, "tau")) %>% 
#   select(m, type, tau_est, tau) %>% 
#   unnest(tau_est, tau) %>% 
#   ggplot(aes(factor(m), tau_est)) +
#   geom_point(colour = 'salmon') +
#   geom_point(aes(y = tau), colour = 'blue') +
#   facet_wrap(~type) + 
#   scale_y_continuous(breaks = 1:15) +
#   labs(y = "Estimated tau", x = "M (number of groups)") +
#   theme_bw()
# 
# ggsave(file = "results/varying_m_tau.png", height = 7, width = 10)


results %>% 
  filter(type == 'stump') %>% 
  mutate(tau_est = map(results_stump, "tau_est"),
         tau = map(results_stump, "tau")) %>% 
  select(m, type, tau_est, tau) %>% 
  unnest(tau_est, tau) %>% 
  ggplot(aes(factor(m), tau_est)) +
  geom_point(colour = 'salmon') +
  geom_point(aes(y = tau), colour = 'blue') +
  facet_wrap(~type) + 
  scale_y_continuous(breaks = 1:15) +
  labs(y = "Estimated tau", x = "M (number of groups)") +
  theme_bw()

ggsave(file = "results/varying_m_stump_tau.png", height = 7, width = 10)
# ------------------------------------------------------------------------------

results %>% 
  filter(type == 'stump') %>% 
  mutate(tau_est = map(results_stump, "tau_est"),
         tau = map(results_stump, "tau"), 
         k1_est = map(results_stump, "k1_est")) %>% 
  select(m, type, tau_est, tau, k1_est) %>% 
  unnest(tau_est, tau, k1_est) %>% 
  mutate(div = tau/k1_est) %>% 
  ggplot(aes(factor(m), div)) +
  geom_point(colour = 'salmon') +
  geom_point(aes(y = tau/8), colour = 'blue') +
  facet_wrap(~type) + 
  labs(y = "Estimated tau/k1", x = "M (number of groups)") +
  theme_bw()

ggsave(file = "results/varying_m_stump_tau_k1.png", height = 7, width = 10)
# ------------------------------------------------------------------------------
# results %>% 
#   mutate(k1_est = map_dbl(results, "k1_est")) %>% 
#   group_by_all() %>% 
#   slice(1) %>% 
#   ggplot(aes(m, k1_est, group = type)) +
#   geom_line(alpha = 0.5) +
#   geom_hline(yintercept = 8) +
#   geom_point(aes(colour = type)) +
#   scale_x_continuous(breaks = c(3, 7, 5, 10, 12, 15, 18, 20, 25, 30)) +
#   scale_y_continuous(breaks = 4:12) +
#   labs(y = "Estimated k1") +
#   theme_bw()
# 
# ggsave(file = "results/varying_m_mean.png", height = 8, width = 7)

# write_rds(results, "results/varying_m.rds")
# 
# m_groups <- 100
# run_100 <- tibble(
#   m = rep(m_groups, 2),
#   type = rep(c("a", "b"), each = length(m_groups))) %>% 
#   mutate(model = map2(m, type, safe_f))

#-------------------------------------------------------------------------------