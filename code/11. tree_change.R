#------------------------------------------------------
#
# September,  2021
#------------------------------------------------------
library(MASS)
library(patchwork)
library(rpart)
library(tidyverse)
source("code/00. simulation_functions.R")
files <- list.files("code/hmb/R") %>% 
  paste0("code/mixedbart/R/", .)

map(files[-5], source)
#------------------------------------------------------
group_variable = "group"
formula <- y ~ X1
iter <-  500
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)


set.seed(2030)
data2 <- sim_b(m = 10, n = 1000)
sqrt(1/data2$tau) 

# rule   <- 0.5
# drawn_node <- "root"
# selec_var <- "X1"
# 
# results_data_b <- data_handler(formula, data2$data, group_variable,
#                                scale_fc = FALSE)
# new_tree <- results_data_b$data %>% 
#   dplyr::mutate(
#     # Increasing the depth of the node giving the grow
#     d =  ifelse(node == drawn_node, d + 1, d),
#     # Updating the parent of the splitted node
#     parent = ifelse(node == drawn_node, drawn_node, parent),
#     # Changing the node "side" of each observation: left and right
#     criteria = ifelse(
#       node == drawn_node,
#       ifelse(!!rlang::sym(selec_var) > rule,
#              "left", "right"), "no split"),
#     
#     # Updating the node accordingly to the new split
#     node = ifelse(node == drawn_node,
#                   ifelse(!!rlang::sym(selec_var) > rule,
#                          paste(node, selec_var, "left"),
#                          paste(node, selec_var, "right")), node),
#     # Updating the node index
#     node_index =  as.numeric(as.factor(node)))

m0_change <- bcart(
  formula, 
  dataset = data2$data,
  iter = 500, 
  group_variable, pars, 
  scale = FALSE, new_tree = TRUE)
saveRDS(m0_change, "results/m0_change.rds")
#current_tree <- m0_change$final_tree

# r <- ratio_grow(tree = new_tree,
#                 current_node =  drawn_node,
#                 pars = pars,
#                 p = 0.5,
#                 current_selec_var = "X1",
#                 results_f = NA,
#                 p_grow = 0.5)


results <- data.frame(
  k1 = m0_change$sampled_k1,
  tau = sqrt(1/m0_change$tau_post[-1])
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

ggsave(file = "results/tau_change.png", height =5, width = 6)

results %>% 
  ggplot(aes(y = k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 8, colour = "red") +
  geom_hline(aes(yintercept = mean(k1)), colour = "blue") +
  labs(y = "Sampled k1s") +
  theme_bw()

ggsave(file = "results/k1_change.png", height = 5, width = 6)

results %>% 
  ggplot(aes(y = tau/k1, iter)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = data2$tau/8, colour = "red") +
  geom_hline(aes(yintercept = mean(tau/k1)), colour = "blue") +
  labs(y = "tau/k1") +
  theme_bw()

ggsave(file = "results/tau_change.png", height = 5, width = 6)

m0_change$final_tree %>% 
  summarise(m = sum((y - mu_sampled)^2))
# 4136
m0_change$final_tree %>% 
  summarise(m = sum((y - mu_js_sampled)^2))

m0_change$final_tree %>%
  group_by(group) %>% 
  summarise(m_y = mean(y), m_mu = mean(mu_js_sampled)) 

m0_change$ratios[200]
m0_change$trees[[201]]$parent

