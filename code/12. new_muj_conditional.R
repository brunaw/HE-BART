#------------------------------------------------------
# This code creates fits a tree by finding the best
# tree and sampling k_1 at the same time 
# August,  2021
#------------------------------------------------------
library(MASS)
library(patchwork)
library(rpart)
library(ranger)
library(tidyverse)
source("code/00. simulation_functions.R")
files <- list.files("code/hmb/R") %>% 
  paste0("code/hmb/R/", .)
map(files, source)
#------------------------------------------------------
group_variable = "group"
formula <- y ~ X1
iter <-  500
#lim <- 50
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)

set.seed(30)
data2 <- sim_b(m = 10, n = 1000)
sqrt(1/data2$tau) #  1.37
sd(data2$data$y)
sd(data2$y2)


m0_regular <- bcart(
  formula, 
  dataset = data2$data,
  iter = iter, 
  group_variable, pars,
  prior_k1 = TRUE,
  min_u = 2, max_u = 12)


n_distinct(m0_regular$final_tree$node)
#saveRDS(m0_regular, "results/m0_full_fit.rds")

results <- data.frame(
  k1 = m0_regular$sampled_k1,
  tau = sqrt(1/m0_regular$tau_post[-1])
) %>% 
  mutate(iter = 1:n())

true_tau <- sqrt(1/data2$tau)

results %>% 
  ggplot(aes(y = tau, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = true_tau, colour = "red") +
  geom_hline(aes(yintercept = mean(tau)), colour = "blue") +
  labs(y = "Sampled taus") +
  theme_bw()

ggsave(file = "results/tau.png", height = 5, width = 6)

results %>% 
  ggplot(aes(y = k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 8, colour = "red") +
  geom_hline(aes(yintercept = mean(k1)), colour = "blue") +
  labs(y = "Sampled k1s") +
  theme_bw()

ggsave(file = "results/k1.png", height = 5, width = 6)
# ----------------------------------------------------------------
n_distinct(m0_regular$final_tree$node)

m0_regular$final_tree %>%
  group_by(group, node) %>% 
  summarise(m_y = mean(y), m_mu = mean(mu_js_sampled))

m0_regular$final_tree %>%
  group_by(node) %>% 
  summarise(m_y = mean(y), m_mu = mean(mu_sampled))

get_rss <- function(model_tree, formula, data,
                    formula2 = y ~ X1 + group){

  rss_hbm <- model_tree %>% 
    summarise(m = sum((y - mu_sampled)^2)) %>% 
    pull(m)
  
  rssj_hbm <- model_tree %>% 
    summarise(m = sum((y - mu_js_sampled)^2)) %>% 
    pull(m)
  
  
  
  b <- rpart::rpart(formula, data, cp = 0.03)
  data$pred_tree <- predict(b)
  
  rss_tree <- data %>% 
    summarise(m = sum((y - pred_tree)^2)) %>% 
    pull(m)
  
  b <- rpart::rpart(formula2, data)
  
  data$pred_tree <- predict(b)
  
  rssj_tree <- data %>% 
    summarise(m = sum((y - pred_tree)^2)) %>% 
    pull(m)
  
  
  b <- ranger(formula, data)
  rss_rf <- sum((b$predictions - data$y)^2)
  b <- ranger(formula2, data)
  # 3287
  rssj_rf <- sum((b$predictions - data$y)^2)
  
  
  df <- data.frame(
  rss = c(rss_hbm, rssj_hbm, rss_tree, rssj_tree, rss_rf, rssj_rf),
  model = rep(c("mhb", "tree", "RF"), each = 2), 
  type = rep(c("RSS without group", "RSS with group"), 3)
)

  
  return(df)
}

df <- get_rss(m0_regular$final_tree, formula = formula, data = data2$data)

df %>%
  spread(type, rss) %>% 
  saveRDS("results/df.rds")
# ----------------------------------------------------------------
iter <-  500
#lim <- 50
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 2, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)

set.seed(55)
data3 <- sim_b(m = 10, n = 1000, k1 = 2)
sqrt(1/data3$tau) #  1.37
sd(data3$data$y)
sd(data3$y2)


m0_regular <- bcart(
  formula, 
  dataset = data3$data,
  iter = iter, 
  group_variable, pars,
  prior_k1 = TRUE,
  min_u = 2, max_u = 12)


n_distinct(m0_regular$final_tree$node)
#saveRDS(m0_regular, "results/m0_full_fit.rds")

results <- data.frame(
  k1 = m0_regular$sampled_k1,
  tau = sqrt(1/m0_regular$tau_post[-1])
) %>% 
  mutate(iter = 1:n())

true_tau <- sqrt(1/data3$tau)
true_tau

results %>% 
  ggplot(aes(y = tau, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = true_tau, colour = "red") +
  geom_hline(aes(yintercept = mean(tau)), colour = "blue") +
  labs(y = "Sampled taus") +
  theme_bw()

ggsave(file = "results/tau_small_k.png", height =5, width = 6)

results %>% 
  ggplot(aes(y = k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, colour = "red") +
  geom_hline(aes(yintercept = mean(k1)), colour = "blue") +
  labs(y = "Sampled k1s") +
  theme_bw()

ggsave(file = "results/k1_small_k.png", height = 5, width = 6)
# ----------------------------------------------------------------
n_distinct(m0_regular$trees[[200]]$node)

m0_regular$final_tree %>%
  group_by(group, node) %>% 
  summarise(m_y = mean(y), m_mu = mean(mu_js_sampled))

m0_regular$final_tree %>%
  group_by(node) %>% 
  summarise(m_y = mean(y), m_mu = mean(mu_sampled))

df <- get_rss(m0_regular$final_tree, 
              formula = formula, data = data3$data)
df %>%
  spread(type, rss) %>% 
  saveRDS("results/df_small_k.rds")
# ----------------------------------------------------------------
# m0_regular <- bcart(
#   formula, 
#   dataset = data2$data,
#   iter = 300, 
#   group_variable, pars,
#   prior_k1 = TRUE,
#   min_u = 2, max_u = 12, 
#   new_tree = TRUE)
# 
# 
# pars$k1 <- 8
# pars$k2 <- 10
# 
# llk <- lk_ratio_grow(
#   m0_regular$final_tree, 
#   current_node = "root",
#   pars
# )
# llk
# map_dbl(m0_regular$trees, ~n_distinct(.x$node))
# current_node <- m0_regular$trees[[4]] %>% 
#   count(node) %>% 
#   slice(1) %>% 
#   pull(node)
# 
# pars$k1 <- 8
# pars$k2 <- 15
# llk <- lk_ratio_grow(
#   tree  = m0_regular$trees[[4]], 
#   current_node = "root X1 left",
#   pars
# )
# llk
# 
# n_distinct(m0_regular$final_tree$node)
# #saveRDS(m0_regular, "results/m0_full_fit.rds")
# 
# results <- data.frame(
#   k1 = m0_regular$sampled_k1,
#   tau = sqrt(1/m0_regular$tau_post[-1])
# ) %>% 
#   mutate(iter = 1:n())
# 
# true_tau <- sqrt(1/data2$tau)
# 
# results %>% 
#   ggplot(aes(y = sqrt(1/tau), iter)) +
#   geom_point(alpha = 0.5) +
#   geom_hline(yintercept = true_tau, colour = "red") +
#   geom_hline(aes(yintercept = mean(sqrt(1/tau))), colour = "blue") +
#   labs(y = "Sampled taus") +
#   theme_bw()
# 
# ggsave(file = "results/tau_true_tree.png", height =5, width = 6)
# 
# results %>% 
#   ggplot(aes(y = k1, iter)) +
#   geom_point(alpha = 0.5) +
#   geom_hline(yintercept = 8, colour = "red") +
#   geom_hline(aes(yintercept = mean(k1)), colour = "blue") +
#   labs(y = "Sampled k1s") +
#   theme_bw()
# 
# ggsave(file = "results/k1_true_tree.png", height = 5, width = 6)
# # ----------------------------------------------------------------
# n_distinct(m0_regular$final_tree$node)
# 
# rss_hbm <- m0_regular$final_tree %>% 
#   summarise(m = sum((y - mu_sampled)^2)) %>% 
#   pull(m)
# rssj_hbm <- m0_regular$final_tree %>% 
#   summarise(m = sum((y - mu_js_sampled)^2)) %>% 
#   pull(m)
# 
# 
# m0_regular$final_tree %>%
#   group_by(group, node) %>% 
#   summarise(m_y = mean(y), m_mu = mean(mu_js_sampled)) 
# 
# b <- rpart::rpart(formula, data2$data,  cp = 0.03)
# rpart.plot::rpart.plot(b)
# 
# data2$data$pred_tree <- predict(b)
# 
# rss_tree <- data2$data %>% 
#   summarise(m = sum((y - pred_tree)^2)) %>% 
#   pull(m)
# 
# b <- rpart::rpart(y ~ X1 + group, data2$data, 
#                   cp = 0.15)
# rpart.plot::rpart.plot(b)
# 
# data2$data$pred_tree <- predict(b)
# 
# rssj_tree <- data2$data %>% 
#   summarise(m = sum((y - pred_tree)^2)) %>% 
#   pull(m)
# 
# library(ranger)s
# b <- ranger(y ~ X1, data2$data)
# rss_rf <- sum((b$predictions - data2$data$y)^2)
# b <- ranger(y ~ X1 + group, data2$data)
# # 3287
# rssj_rf <- sum((b$predictions - data2$data$y)^2)
# 
# 
# df <- data.frame(
#   rss = c(rss_hbm, rssj_hbm, rss_tree, rssj_tree, rss_rf, rssj_rf),
#   model = rep(c("mhb", "tree", "RF"), each = 2), 
#   type = rep(c("RSS without group", "RSS with group"), 3)
# )
# 
# 
# df %>% 
#   spread(type, rss) 
# #write.table("results/rss.txt")
