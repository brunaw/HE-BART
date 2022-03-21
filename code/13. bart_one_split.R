#-------------------------------------------------------------------------------
# A Hierarchical BART model run for a one-split, fixed tree  
# 


#-------------------------------------------------------------------------------
library(tidyverse)
library(ranger)
source("code/00. simulation_functions.R")

files <- list.files("code/hmb/R") %>% 
  paste0("code/hmb/R/", .)
map(files, source)
#-------------------------------------------------------------------------------
set.seed(100)
data_sim <- sim_bart(n = 1000, m = 10, k2 = 5)
data_sim$tau # tau = 0.325548
group_variable = "group"
formula <- y ~ X1
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 5, alpha = alpha, beta = beta, mu_mu = 0
)

# -------------------------------------------------------------------------
# do we simulate the same way -- no 
m0_2node <- bart(formula, 
                 dataset = data_sim$data,
                 iter = 150, 
                 group_variable, pars, 
                 scale = FALSE)

df_tau <- data.frame(tau = m0_2node$tau_post) %>% 
  slice(-1)

tau_TRUE <- data_sim$tau

p22 <- ggplot(df_tau, aes(
  x = tau
)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = tau_TRUE, colour = 'red') +
  ggtitle("tau posterior") +
  theme_bw()
p22
ggsave(file = "results/reports/december/tau_fixed_one_node.png")

cc <- m0_2node$final_trees %>% 
  mutate(current_tree = map(tree_data, ~{ slice(.x, (150))})) %>% 
  select(tree_index, current_tree) %>% 
  unnest(current_tree) %>% 
  select(tree_index, est_tree) %>% 
  unnest(est_tree) %>% 
  dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>% 
  ungroup()


res <- cc %>% 
  dplyr::select(tree_index, y, sampled_mu_j) %>% 
  pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>% 
  mutate(pred = rowSums(.[2:ncol(.)]), 
         res = (y - pred)**2) 

res_mu <- cc %>% 
  dplyr::select(tree_index, y, sampled_mu) %>% 
  pivot_wider(names_from = tree_index, values_from = sampled_mu) %>% 
  mutate(pred = rowSums(.[2:ncol(.)]), 
         res = (y - pred)**2) 


res %>% 
  ggplot(aes(y, pred)) +
  geom_point() +
  labs(x = 'y', y = 'prediction') +
  theme_bw()

ggsave(file = "results/reports/december/predictions_fixed_one_node.png")

# res_muj <- cc %>% 
#   group_by(tree_index, group, node) %>% 
#   summarise(sampled_mu_j = unique(sampled_mu_j),
#             sampled_mu = unique(sampled_mu), 
#             res = (sampled_mu_j - sampled_mu)**2)
# 
# res_mu <- cc %>% 
#   group_by(tree_index, node) %>% 
#   summarise(sampled_mu = first(sampled_mu), 
#             res = (sampled_mu)**2)

rss <- function(x, y) sum((x - y)**2)
rss_df <- data.frame(rss = rep(NA, 6), source = rep(NA, 6))
# rpart
m0 <- rpart::rpart(y ~ X1, data_sim$data)
m1 <- rpart::rpart(y ~ X1 + group, data_sim$data)
rss_df$rss[1] <- rss(data_sim$data$y, predict(m0, data_sim$data))
rss_df$rss[2]  <- rss(data_sim$data$y, predict(m1, data_sim$data))

# rf
m0 <- ranger::ranger(y ~ X1, data_sim$data, num.trees = 5)
m1 <- ranger::ranger(y ~ X1 + group, data_sim$data, num.trees = 5)
rss_df$rss[3] <- rss(data_sim$data$y, predict(m0,  data_sim$data)$predictions)
rss_df$rss[4] <- rss(data_sim$data$y, predict(m1,  data_sim$data)$predictions)
rss_df$rss[6] <- sum(res$res)


rss_df$rss[5] <- sum(res_mu$res)

rss_df$source <- c("tree", "tree", "RF", "RF", "HBart", "HBart")
rss_df$type <- rep(c("normal", "group"), 3)

rss_df %>% 
  select(source, type, rss) %>% 
  write.table(file = "results/reports/december/rss.txt")

# > m0 <- rpart::rpart(y ~ X1, data_sim$data)
# > m1 <- rpart::rpart(y ~ X1 + group, data_sim$data)
# > rss(data_sim$data$y, predict(m0, data_sim$data))
# [1] 9201.931
# > rss(data_sim$data$y, predict(m1, data_sim$data))
# [1] 1209.953
# > # rf
#   > m0 <- ranger::ranger(y ~ X1, data_sim$data, num.trees = 5)
# > m1 <- ranger::ranger(y ~ X1 + group, data_sim$data, num.trees = 5)
# > rss(data_sim$data$y, predict(m0,  data_sim$data)$predictions)
# [1] 3925.085
# > rss(data_sim$data$y, predict(m1,  data_sim$data)$predictions)
# [1] 1225.62
# > sum(res$res)
# [1] 1175.315
