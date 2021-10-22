#------------------------------------------------------
# This code creates fits a tree by finding the best
# tree and sampling k_1 at the same time 
# August,  2021
#------------------------------------------------------
library(MASS)
library(patchwork)
library(rpart)
library(tidyverse)
source("code/00. simulation_functions.R")
files <- list.files("code/mixedbart/R") %>% 
  paste0("code/mixedbart/R/", .)
map(c(files[-c(7, 10, 11)]), source)
#------------------------------------------------------
group_variable = "group"
formula <- y ~ X1
iter <-  2500
lim <- 50
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)


set.seed(2030)
data2 <- sim_b(m = 10, n = 2000)
sqrt(1/data2$tau) #  1.230742

m0_regular <- bcart(
  formula, 
  dataset = data2$data,
  iter = iter, 
  group_variable, pars)


n_distinct(m0_regular$final_tree$node)
#saveRDS(m0_regular, "results/m0_full_fit.rds")

results <- data.frame(
  k1 = m0_regular$sampled_k1,
  tau = sqrt(1/m0_regular$tau_post[-1])
) %>% 
  mutate(iter = 1:n())

true_tau <- sqrt(1/data2$tau)

results %>% 
  ggplot(aes(y = sqrt(1/tau), iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = true_tau, colour = "red") +
  geom_hline(aes(yintercept = mean(sqrt(1/tau))), colour = "blue") +
  labs(y = "Sampled taus") +
  theme_bw()

#ggsave(file = "results/tau_full_fit.png", height =5, width = 6)

results %>% 
  ggplot(aes(y = k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 8, colour = "red") +
  geom_hline(aes(yintercept = mean(k1)), colour = "blue") +
  labs(y = "Sampled k1s") +
  theme_bw()

#ggsave(file = "results/k1_full_fit.png", height = 5, width = 6)

results %>% 
  ggplot(aes(y = tau/k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = data2$tau/8, colour = "red") +
  geom_hline(aes(yintercept = mean(tau/k1)), colour = "blue") +
  labs(y = "tau/k1") +
  theme_bw()


#ggsave(file = "results/tau_k1_full_fit.png", height = 5, width = 6)
# ----------------------------------------------------------------
n_distinct(m0_regular$final_tree$node)

# 26200
m0_regular$final_tree %>% 
  summarise(m = sum((y - mu_sampled)^2))
# 4136
m0_regular$final_tree %>% 
  summarise(m = sum((y - mu_js_sampled)^2))

m0_regular$final_tree %>%
  group_by(group) %>% 
  summarise(m_y = mean(y), m_mu = mean(mu_js_sampled)) 

b <- rpart::rpart(formula, data2$data,  cp = 0.03)
rpart.plot::rpart.plot(b)

m0_regular$final_tree$pred_tree <- predict(b)

m0_regular$final_tree %>% 
  summarise(m = sum((y - pred_tree)^2))

sum((predict(b) - data2$data$y)^2)

b <- rpart::rpart(y ~ X1 + group, data2$data, 
                  cp = 0.15)
rpart.plot::rpart.plot(b)

m0_regular$final_tree$pred_tree <- predict(b)

m0_regular$final_tree %>% 
  summarise(m = sum((y - pred_tree)^2))

m0_regular$final_tree %>% 
  group_by(group) %>% 
  summarise(m_y = mean(y), m_mu = mean(pred_tree))

m0_regular$final_tree %>% 
  ggplot(aes(x = y)) +
  geom_point(aes(y = mu_js_sampled)) +
  #geom_point(aes(y = pred_tree), colour = "red") +
  facet_wrap(~group) +
  theme_bw()

library(ranger)
b <- ranger(y ~ X1, data2$data)
sum((b$predictions - data2$data$y)^2)
b <- ranger(y ~ X1 + group, data2$data)
# 3287
sum((b$predictions - data2$data$y)^2)
# ----------------------------------------------------------------
# running more iterations
iter <-  5000
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)

m0_regular_iter <- bcart(
  formula, 
  dataset = data2$data,
  iter = iter, 
  group_variable, pars)

n_distinct(m0_regular_iter$final_tree$node)

results <- data.frame(
  k1 = m0_regular_iter$sampled_k1,
  tau = sqrt(1/m0_regular_iter$tau_post[-1])
) %>% 
  mutate(iter = 1:n())

true_tau <- sqrt(1/data2$tau)

results %>% 
  ggplot(aes(y = sqrt(1/tau), iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = true_tau, colour = "red") +
  geom_hline(aes(yintercept = mean(sqrt(1/tau))), colour = "blue") +
  labs(y = "Sampled taus") +
  theme_bw()

ggsave(file = "results/tau_full_fit_iter.png", 
       height =5, width = 6)

results %>% 
  ggplot(aes(y = k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 8, colour = "red") +
  geom_hline(aes(yintercept = mean(k1)), colour = "blue") +
  labs(y = "Sampled k1s") +
  theme_bw()

ggsave(file = "results/k1_full_fit_iter.png", height = 5, width = 6)

results %>% 
  ggplot(aes(y = tau/k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = true_tau/8, colour = "red") +
  geom_hline(aes(yintercept = mean(tau/k1)), colour = "blue") +
  labs(y = "tau/k1") +
  theme_bw()

ggsave(file = "results/tau_k1_full_fit_iter.png", height = 5, width = 6)

saveRDS(m0_regular_iter, "results/m0_full_fit_iter.rds")
# ----------------------------------------------------------------
m0_regular_iter$final_tree %>% 
  summarise(m = sum((y - mu_sampled)^2))

m0_regular_iter$final_tree %>% 
  summarise(m = sum((y - mu_js_sampled)^2))

# m0_regular_iter$final_tree %>%
#   group_by(node, group) %>% 
#   summarise(m_y = mean(y), m_mu = mean(mu_js_sampled)) %>% 
#   View()

m0_regular$final_tree %>%
  group_by(group) %>% 
  summarise(m_y = mean(y), m_mu = mean(mu_sampled)) 

b <- rpart::rpart(formula, data2$data)
rpart.plot::rpart.plot(b)
data2$data$pred_tree <- predict(b)

data2$data %>% 
  group_by(group) %>% 
  summarise(m_y = mean(y), m_mu = mean(pred_tree))

sum((predict(b) - data2$data$y)^2)

b <- rpart::rpart(y ~ X1 + group, data2$data)
sum((predict(b) - data2$data$y)^2)

# ----------------------------------------------------------------
# changing the value of k1
iter = 2500
set.seed(2021)
data3 <- sim_b(m = 5, n = 2000, k1 = 12)
sqrt(1/data3$tau) 

m0_regular_12 <- bcart(
  formula, 
  dataset = data3$data,
  iter = iter, 
  group_variable, pars)

n_distinct(m0_regular_12$final_tree$node)
saveRDS(m0_regular_12, "results/m0_full_fit_12.rds")
# ----------------------------------------------------------------
results <- data.frame(
  k1 = m0_regular_12$sampled_k1,
  tau = sqrt(1/m0_regular_12$tau_post[-1])
) %>% 
  mutate(iter = 1:n())

true_tau <- sqrt(1/data3$tau)

results %>% 
  ggplot(aes(y = sqrt(1/tau), iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = true_tau, colour = "red") +
  geom_hline(aes(yintercept = mean(sqrt(1/tau))), colour = "blue") +
  labs(y = "Sampled taus") +
  theme_bw()

ggsave(file = "results/tau_full_fit_12.png", height =5, width = 6)

results %>% 
  ggplot(aes(y = k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 12, colour = "red") +
  geom_hline(aes(yintercept = mean(k1)), colour = "blue") +
  labs(y = "Sampled k1s") +
  theme_bw()

ggsave(file = "results/k1_full_fit_12.png", height = 5, width = 6)

results %>% 
  ggplot(aes(y = tau/k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = data3$tau/12, colour = "red") +
  geom_hline(aes(yintercept = mean(tau/k1)), colour = "blue") +
  labs(y = "tau/k1") +
  theme_bw()


ggsave(file = "results/tau_k1_full_fit_12.png", height = 5, width = 6)
# ----------------------------------------------------------------
png(filename = "results/prior_k1.png")
curve(dweibull(x, 7, 10), 0, 15, col = '2')
dev.off()

# Adding the prior 
group_variable = "group"
formula <- y ~ X1
iter <-  2500
lim <- 200
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 5, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)

set.seed(2025)
data2 <- sim_b(m = 10, n = 2000)
sqrt(1/data2$tau) 

m0_regular_prior <- bcart(
  formula, 
  dataset = data2$data,
  iter = iter, 
  group_variable, pars, 
  prior_k1 = TRUE,
  min_u = 2, max_u = 12)

saveRDS(m0_regular_prior, "results/m0_full_fit_prior.rds")
#-----------------------------------------------------------------------
# Prior evaluation
current_tree <- m0_regular_prior$final_tree
n_distinct(current_tree$node)

k1_seq <- seq(0.01, 25, by = 0.1)
ll_k1 <- map_dbl(k1_seq, calc_all, current_tree = current_tree, pars = pars)
ll_prior <- dweibull(k1_seq, shape = 7, 10, log = TRUE)

p1 <- data.frame(k1_seq, ll_k1) %>% 
  ggplot(aes(k1_seq, ll_k1)) +
  geom_point() +
  labs(y = "Conditional", x = "k1 values",
       title = "Results using stump") +
  theme_bw() 


p2 <- data.frame(k1_seq, ll_prior) %>% 
  ggplot(aes(k1_seq, ll_prior)) +
  geom_point() +
  labs(y = "Prior (dweibull(7, 10))", x = "k1 values") +
  theme_bw() 


p3 <- data.frame(k1_seq, ll = ll_k1, ll_prior = ll_prior) %>% 
  mutate(ll_2 = ll - lag(ll) + ll_prior) %>% 
  ggplot(aes(k1_seq, ll_2)) +
  geom_point() +
  labs(y = "Ratio (ll - lag(ll)) + prior", x = "k1 values") +
  theme_bw() 


p4 <- data.frame(k1_seq, ll = ll_k1, ll_prior = ll_prior) %>% 
  mutate(ll_2 = ll - lag(ll) + ll_prior, 
         exp = exp(ll_2)) %>% 
  filter(k1_seq > 1.2) %>% 
  ggplot(aes(k1_seq, exp)) +
  geom_point() +
  labs(y = "Exp(ll ratio + prior)", 
       x = "k1 values") +
  theme_bw() 

p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
ggsave(file = "results/ll_exp_fit_prior.png", height = 6, width = 8)
#-----------------------------------------------------------------------
results <- data.frame(
  k1 = m0_regular_prior$sampled_k1,
  tau = sqrt(1/m0_regular_prior$tau_post[-1])
) %>% 
  mutate(iter = 1:n())

true_tau <- sqrt(1/data2$tau)

results %>% 
  ggplot(aes(y = sqrt(1/tau), iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = true_tau, colour = "red") +
  geom_hline(aes(yintercept = mean(sqrt(1/tau))), colour = "blue") +
  labs(y = "Sampled taus") +
  theme_bw()

ggsave(file = "results/tau_full_fit_prior.png", height =5, width = 6)

results %>% 
  ggplot(aes(y = k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 8, colour = "red") +
  geom_hline(aes(yintercept = mean(k1)), colour = "blue") +
  labs(y = "Sampled k1s") +
  theme_bw()

ggsave(file = "results/k1_full_fit_prior.png", height = 5, width = 6)

results %>% 
  ggplot(aes(y = tau/k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = data2$tau/8, colour = "red") +
  geom_hline(aes(yintercept = mean(tau/k1)), colour = "blue") +
  labs(y = "tau/k1") +
  theme_bw()


ggsave(file = "results/tau_k1_full_fit_prior.png", height = 5, width = 6)
# ----------------------------------------------------------------
# Stopped here
# ----------------------------------------------------------------
