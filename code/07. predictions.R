#-----------------------------------------------------------------
all_data <- readRDS("data/all_data.rds")
data_one_node <- readRDS("data/data_one_node.rds")
data_k1 <- readRDS("data/data_k1.rds")
library(tidyverse)
library(rpart)
library(ranger)
files <- list.files("code/mixedbart/R") %>% 
  paste0("code/mixedbart/R/", .)
map(c(files[-c(7, 10, 11)]), source)

rss <- function(y, y_test){ mean((y - y_test)^2)}
#-----------------------------------------------------------------
m0_2node          <- readRDS("results/m0_2node.rds")
m0_3node          <- readRDS("results/m0_3node.rds")
m1_3node          <- readRDS("results/m1_3node.rds")
m2_3node          <- readRDS("results/m2_3node.rds")
m3_3node          <- readRDS("results/m3_3node.rds")
m0_regular        <- readRDS("results/m0_regular.rds")
m0_k1_comp        <- readRDS("results/m0_k1_comp.rds")
m0_k1_samp        <- readRDS("results/m0_k1_samp.rds")
m0_k1_samp_bigger <- readRDS("results/m0_k1_bigger.rds")
#-----------------------------------------------------------------
# rmse_results <- data.frame(
#   model = c("2 node", "3 node, small ks", 
#             "3 node, big k2", "3 node, big k1", 
#             "3 node, big ks", "3 node, non-fixed tree"))
# rmse_results$hbcart <- "NA"
# rmse_results$hbcart_group <- "NA"
# rmse_results$tree <- "NA"
# rmse_results$tree_group <- "NA"
# rmse_results$rf <- "NA"
# rmse_results$rf_group <- "NA"

rmse_results <- read.table("results/rmse_results.txt")
#----------------------------------------------------------------- 
final_mus <- res_2node %>% 
  select(2:3) %>% 
  gather() %>% 
  mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  summarise(mean_sample = mean(value)) %>% 
  mutate(node = unique(m0_2node$final_tree$node)) %>% 
  select(-key)

final_means <- m0_2node$final_tree %>%
  group_by(node, group) %>% 
  summarise(final_mu = mean(mu_js_sampled), 
            y = mean(y)) %>% 
  left_join(final_mus, by = "node")

rmse_one_node <- data_one_node$test %>% 
  mutate(
    node = 
      ifelse(X1 > 0.5,"root X1 left", "root X1 right")) %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  summarise(rmse = mean((y - final_mu)^2), 
            rmse_mu = mean((y - mean_sample)^2))  


rmse_results[1, 2] <- rmse_one_node$rmse_mu
rmse_results[1, 3] <- rmse_one_node$rmse

data_one_node$test %>% 
  mutate(
    node = 
      ifelse(X1 > 0.5,"root X1 left", "root X1 right")) %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  group_by(node, group) %>% 
  summarise(final_mu = mean(final_mu), 
            y = mean(y))

# rpart
m0 <- rpart::rpart(y ~ X1, data_one_node$train)
m1 <- rpart::rpart(y ~ X1 + group, data_one_node$train)
rss(data_one_node$train$y, predict(m0, data_one_node$train))
rss(data_one_node$train$y, predict(m1, data_one_node$train))

rmse_results[1, 4] <- rss(data_one_node$test$y, predict(m0, data_one_node$test))
rmse_results[1, 5] <- rss(data_one_node$test$y, predict(m1, data_one_node$test))

# rf
m0 <- ranger::ranger(y ~ X1, data_one_node$train, num.trees = 100)
m1 <- ranger::ranger(y ~ X1 + group, data_one_node$train, num.trees = 100)
rmse_results[1, 6] <- rss(data_one_node$test$y, predict(m0,  data_one_node$test)$predictions)
rmse_results[1, 7] <- rss(data_one_node$test$y, predict(m1,  data_one_node$test)$predictions)
#---------------------------------------------------------------------
final_mus <- res_3node %>% 
  select(2:4) %>% 
  gather() %>% 
  mutate(true_mean = rep(sim_ks_small$mu[c(1,2,3)], 
                         each = (iter-lim+1))) %>%
  group_by(key) %>% 
  summarise(mean_sample = mean(value)) %>% 
  select(-key) %>% 
  mutate(node = unique(m0_3node$final_tree$node))

final_means <- m0_3node$final_tree %>%
  group_by(node, group) %>% 
  summarise(final_mu = mean(mu_js_sampled), 
            y = mean(y)) %>% 
  left_join(final_mus, by = "node")

# final_means %>%
#   group_by(node) %>%
#   summarise(m = mean(final_mu), m1 = mean(mean_sample),
#             y = mean(y))

df <- sim_ks_small$data$test %>% 
  mutate(
    node = ifelse(X1 > 0.5, "root X1 left",
                  "root X1 right X2 left"), 
    node = ifelse(
      node == "root X1 right X2 left" & X2 < 0.5, 
      "root X1 right X2 right", node)) 

rmse <- df %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  summarise(rmse = mean((y - final_mu)^2), 
            rmse_mu = mean((y - mean_sample)^2))  


rmse_results[2, 2] <- rmse$rmse_mu
rmse_results[2, 3] <- rmse$rmse

# rpart
m0 <- rpart::rpart(y ~ X1 + X2, sim_ks_small$data$train)
m1 <- rpart::rpart(y ~ X1 + X2 + group, sim_ks_small$data$train)
rmse_results[2, 4] <- rss(sim_ks_small$data$test$y, predict(m0, sim_ks_small$data$test))
rmse_results[2, 5] <- rss(sim_ks_small$data$test$y, predict(m1, sim_ks_small$data$test))

# rf
m0 <- ranger::ranger(y ~ X1 + X2, sim_ks_small$data$train, num.trees = 100)
m1 <- ranger::ranger(y ~ X1 + X2 + group, sim_ks_small$data$train, num.trees = 100)
rmse_results[2, 6] <- rss(sim_ks_small$data$test$y, 
    predict(m0, sim_ks_small$data$test)$predictions)
rmse_results[2, 7] <- rss(sim_ks_small$data$test$y, 
    predict(m1, sim_ks_small$data$test)$predictions)
#---------------------------------------------------------------------
final_mus <- res_3node2 %>% 
  select(2:4) %>% 
  gather() %>% 
  mutate(true_mean = rep(sim_k2_big$mu[c(1,2,3)], 
                         each = (iter-lim+1))) %>%
  group_by(key) %>% 
  summarise(mean_sample = mean(value)) %>% 
  select(-key) %>% 
  mutate(node = unique(m1_3node$final_tree$node))

final_means <- m1_3node$final_tree %>%
  group_by(node, group) %>% 
  summarise(final_mu = mean(mu_js_sampled), 
            y = mean(y)) %>% 
  left_join(final_mus, by = "node")

# final_means %>%
#   group_by(node) %>%
#   summarise(m = mean(final_mu), m1 = mean(mean_sample),
#             y = mean(y))

df <- sim_k2_big$data$test %>% 
  mutate(
    node = ifelse(X1 > 0.5, "root X1 left",
                  "root X1 right X2 left"), 
    node = ifelse(
      node == "root X1 right X2 left" & X2 < 0.5, 
      "root X1 right X2 right", node)) 

rmse <- df %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  summarise(rmse = mean((y - final_mu)^2), 
            rmse_mu = mean((y - mean_sample)^2))  


rmse_results[3, 2] <- rmse$rmse_mu
rmse_results[3, 3] <- rmse$rmse


# rpart
m0 <- rpart::rpart(y ~ X1 + X2, sim_k2_big$data$train)
m1 <- rpart::rpart(y ~ X1 + X2 + group, sim_k2_big$data$train)
rmse_results[3, 4] <- rss(sim_k2_big$data$test$y, predict(m0, sim_k2_big$data$test))
rmse_results[3, 5] <- rss(sim_k2_big$data$test$y, predict(m1, sim_k2_big$data$test))

# rf
m0 <- ranger::ranger(y ~ X1 + X2, sim_k2_big$data$train, num.trees = 100)
m1 <- ranger::ranger(y ~ X1 + X2 + group, sim_k2_big$data$train, num.trees = 100)
rmse_results[3, 6] <- rss(sim_k2_big$data$test$y, predict(m0,  
                                    sim_k2_big$data$test)$predictions)
rmse_results[3, 7] <- rss(sim_k2_big$data$test$y, predict(m1,  
                                    sim_k2_big$data$test)$predictions)
#---------------------------------------------------------------------
final_mus <- res_3node3 %>% 
  select(2:4) %>% 
  gather() %>% 
  mutate(true_mean = rep(sim_k1_big$mu[c(2,1,3)], 
                         each = (iter-lim+1))) %>%
  group_by(key) %>% 
  summarise(mean_sample = mean(value)) %>% 
  select(-key) %>% 
  mutate(node = unique(m2_3node$final_tree$node))

final_means <- m2_3node$final_tree %>%
  group_by(node, group) %>% 
  summarise(final_mu = mean(mu_js_sampled), 
            y = mean(y)) %>% 
  left_join(final_mus, by = "node")

# final_means %>%
#   group_by(node) %>%
#   summarise(m = mean(final_mu), m1 = mean(mean_sample),
#             y = mean(y))

df <- sim_k1_big$data$test %>% 
  mutate(
    node = ifelse(X1 > 0.5, "root X1 left",
                  "root X1 right X2 left"), 
    node = ifelse(
      node == "root X1 right X2 left" & X2 < 0.5, 
      "root X1 right X2 right", node)) 

rmse <- df %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  summarise(rmse = mean((y - final_mu)^2), 
            rmse_mu = mean((y - mean_sample)^2))  


rmse_results[4, 2] <- rmse$rmse_mu
rmse_results[4, 3] <- rmse$rmse

# rpart
m0 <- rpart::rpart(y ~ X1 + X2, sim_k1_big$data$train)
m1 <- rpart::rpart(y ~ X1 + X2 + group, sim_k1_big$data$train)
rmse_results[4, 4] <- rss(sim_k1_big$data$test$y, predict(m0, sim_k1_big$data$test))
rmse_results[4, 5] <- rss(sim_k1_big$data$test$y, predict(m1, sim_k1_big$data$test))

# rf
m0 <- ranger::ranger(y ~ X1 + X2, sim_k1_big$data$train, num.trees = 100)
m1 <- ranger::ranger(y ~ X1 + X2 + group, sim_k1_big$data$train, num.trees = 100)
rmse_results[4, 6] <- rss(sim_k1_big$data$test$y, predict(m0,  
                                    sim_k1_big$data$test)$predictions)
rmse_results[4, 7] <- rss(sim_k1_big$data$test$y, predict(m1,  
                                    sim_k1_big$data$test)$predictions)
#---------------------------------------------------------------------
final_mus <- res_3node4 %>% 
  select(2:4) %>% 
  gather() %>% 
  mutate(true_mean = rep(sim_ks_big$mu[c(1,2,3)], 
                         each = (iter-lim+1))) %>%
  group_by(key) %>% 
  summarise(mean_sample = mean(value)) %>% 
  select(-key) %>% 
  mutate(node = unique(m3_3node$final_tree$node))

final_means <- m3_3node$final_tree %>%
  group_by(node, group) %>% 
  summarise(final_mu = mean(mu_js_sampled), 
            y = mean(y)) %>% 
  left_join(final_mus, by = "node")

# final_means %>%
#   group_by(node) %>%
#   summarise(m = mean(final_mu), m1 = mean(mean_sample),
#             y = mean(y))

df <- sim_ks_big$data$test %>% 
  mutate(
    node = ifelse(X1 > 0.5, "root X1 left",
                  "root X1 right X2 left"), 
    node = ifelse(
      node == "root X1 right X2 left" & X2 < 0.5, 
      "root X1 right X2 right", node)) 

rmse <- df %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  summarise(rmse = mean((y - final_mu)^2), 
            rmse_mu = mean((y - mean_sample)^2))  


rmse_results[5, 2] <- rmse$rmse_mu
rmse_results[5, 3] <- rmse$rmse

# rpart
m0 <- rpart::rpart(y ~ X1 + X2, sim_ks_big$data$train)
m1 <- rpart::rpart(y ~ X1 + X2 + group, sim_ks_big$data$train)
rmse_results[5, 4] <- rss(sim_ks_big$data$test$y, predict(m0, sim_ks_big$data$test))
rmse_results[5, 5] <- rss(sim_ks_big$data$test$y, predict(m1, sim_ks_big$data$test))

# rf
m0 <- ranger::ranger(y ~ X1 + X2, sim_ks_big$data$train, num.trees = 100)
m1 <- ranger::ranger(y ~ X1 + X2 + group, sim_ks_big$data$train, num.trees = 100)
rmse_results[5, 6] <- rss(sim_ks_big$data$test$y, predict(m0,  
                                    sim_ks_big$data$test)$predictions)
rmse_results[5, 7] <- rss(sim_ks_big$data$test$y, predict(m1,  
                                    sim_ks_big$data$test)$predictions)

#---------------------------------------------------------------------
final_mus <- res_regular %>% 
  select(2:3) %>% 
  gather() %>% 
  mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  summarise(mean_sample = mean(value)) %>% 
  mutate(node = unique(m0_regular$final_tree$node)) %>% 
  select(-key)

final_means <- m0_regular$final_tree %>%
  group_by(node, group) %>% 
  summarise(final_mu = mean(mu_js_sampled), 
            y = mean(y)) %>% 
  left_join(final_mus, by = "node")

# final_means %>%
#   group_by(node) %>%
#   summarise(m = mean(final_mu), m1 = mean(mean_sample),
#             y = mean(y))

df <- sim_ks_big$data$test %>% 
  mutate(
    node = ifelse(X1 < 0.595, "root X1 right",
                  "root X1 left")) 


rmse <- df %>%                                 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) 


rmse <- rmse %>% 
  summarise(rmse = mean((y - final_mu)^2), 
            rmse_mu = mean((y - mean_sample)^2))  

rmse_results[6, 2] <- rmse$rmse_mu
rmse_results[6, 3] <- rmse$rmse


# rpart
m0 <- rpart::rpart(y ~ X1 + X2, sim_ks_big$data$train)
m1 <- rpart::rpart(y ~ X1 + X2 + group, sim_ks_big$data$train)
rmse_results[6, 4] <- rss(sim_ks_big$data$test$y, predict(m0, sim_ks_big$data$test))
rmse_results[6, 5] <- rss(sim_ks_big$data$test$y, predict(m1, sim_ks_big$data$test))

# rf
m0 <- ranger::ranger(y ~ X1 + X2, sim_ks_big$data$train, num.trees = 100)
m1 <- ranger::ranger(y ~ X1 + X2 + group, sim_ks_big$data$train, num.trees = 100)
rmse_results[6, 6] <- rss(sim_ks_big$data$test$y, predict(m0,  
                                    sim_ks_big$data$test)$predictions)
rmse_results[6, 7] <- rss(sim_ks_big$data$test$y, predict(m1,  
                                    sim_ks_big$data$test)$predictions)

#-----------------------------------------------------------------
mu = c(-1.628563, 1.483883)

final_mus <- res_sample_k1 %>% 
  select(2:3) %>% 
  gather() %>% 
  mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  summarise(mean_sample = mean(value)) %>% 
  mutate(node = unique(m0_k1_samp$final_tree$node)) %>% 
  select(-key)

final_means <- m0_k1_samp$final_tree %>%
  group_by(node, group) %>% 
  summarise(final_mu = mean(mu_js_sampled), 
            y = mean(y)) %>% 
  left_join(final_mus, by = "node")

rmse_k1 <- data_k1$test %>% 
  mutate(
    node = 
      ifelse(X1 > 0.5,"root X1 left", "root X1 right")) %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  summarise(rmse = mean((y - final_mu)^2), 
            rmse_mu = mean((y - mean_sample)^2))  


rmse_results[7, 2] <- rmse_k1$rmse_mu
rmse_results[7, 3] <- rmse_k1$rmse

data_k1$test %>% 
  mutate(
    node = 
      ifelse(X1 > 0.5,"root X1 left", "root X1 right")) %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  group_by(node, group) %>% 
  summarise(final_mu = mean(final_mu), 
            y = mean(y))

# rpart
m0 <- rpart::rpart(y ~ X1, data_k1$train)
m1 <- rpart::rpart(y ~ X1 + group, data_k1$train)
rss(data_k1$train$y, predict(m0, data_k1$train))
rss(data_k1$train$y, predict(m1, data_k1$train))

rmse_results[7, 4] <- rss(data_k1$test$y, predict(m0, data_k1$test))
rmse_results[7, 5] <- rss(data_k1$test$y, predict(m1, data_k1$test))

# rf
m0 <- ranger::ranger(y ~ X1, data_k1$train, num.trees = 100)
m1 <- ranger::ranger(y ~ X1 + group, data_k1$train, num.trees = 100)
rmse_results[7, 6] <- rss(data_k1$test$y, predict(m0,  data_k1$test)$predictions)
rmse_results[7, 7] <- rss(data_k1$test$y, predict(m1,  data_k1$test)$predictions)
#rmse_results$model <- as.character(rmse_results$model)
rmse_results[7, 1] <- "2 node, sampling k1"
#-----------------------------------------------------------------------
final_mus <- res_sample_comp %>% 
  select(2:3) %>% 
  gather() %>% 
  mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  summarise(mean_sample = mean(value)) %>% 
  mutate(node = unique(m0_k1_comp$final_tree$node)) %>% 
  select(-key)

final_means <- m0_k1_comp$final_tree %>%
  group_by(node, group) %>% 
  summarise(final_mu = mean(mu_js_sampled), 
            y = mean(y)) %>% 
  left_join(final_mus, by = "node")

rmse_comp <- data_k1$test %>% 
  mutate(
    node = 
      ifelse(X1 > 0.5,"root X1 left", "root X1 right")) %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  summarise(rmse = mean((y - final_mu)^2), 
            rmse_mu = mean((y - mean_sample)^2))  


rmse_results[8, 2] <- rmse_comp$rmse_mu
rmse_results[8, 3] <- rmse_comp$rmse

data_k1$test %>% 
  mutate(
    node = 
      ifelse(X1 > 0.5,"root X1 left", "root X1 right")) %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  group_by(node, group) %>% 
  summarise(final_mu = mean(final_mu), 
            y = mean(y))

# rpart
m0 <- rpart::rpart(y ~ X1, data_k1$train)
m1 <- rpart::rpart(y ~ X1 + group, data_k1$train)
rss(data_k1$train$y, predict(m0, data_k1$train))
rss(data_k1$train$y, predict(m1, data_k1$train))

rmse_results[8, 4] <- rss(data_k1$test$y, predict(m0, data_k1$test))
rmse_results[8, 5] <- rss(data_k1$test$y, predict(m1, data_k1$test))

# rf
m0 <- ranger::ranger(y ~ X1, data_k1$train, num.trees = 100)
m1 <- ranger::ranger(y ~ X1 + group, data_k1$train, num.trees = 100)
rmse_results[8, 6] <- rss(data_k1$test$y, predict(m0,  data_k1$test)$predictions)
rmse_results[8, 7] <- rss(data_k1$test$y, predict(m1,  data_k1$test)$predictions)
rmse_results$model <- as.character(rmse_results$model)
rmse_results$model[8] <- "2 node, same data as previous one"
#-----------------------------------------------------------------
final_mus <- res_sample_k1_big %>% 
  select(2:3) %>% 
  gather() %>% 
  mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  summarise(mean_sample = mean(value)) %>% 
  mutate(node = unique(m0_k1_samp_bigger$final_tree$node)) %>% 
  select(-key)

final_means <- m0_k1_samp_bigger$final_tree %>%
  group_by(node, group) %>% 
  summarise(final_mu = mean(mu_js_sampled), 
            y = mean(y)) %>% 
  left_join(final_mus, by = "node")

rmse_k1_big <- data_k1$test %>% 
  mutate(
    node = 
      ifelse(X1 > 0.5,"root X1 left", "root X1 right")) %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  summarise(rmse = mean((y - final_mu)^2), 
            rmse_mu = mean((y - mean_sample)^2))  


rmse_results[9, 2] <- rmse_k1_big$rmse_mu
rmse_results[9, 3] <- rmse_k1_big$rmse

data_k1$test %>% 
  mutate(
    node = 
      ifelse(X1 > 0.5,"root X1 left", "root X1 right")) %>% 
  left_join(final_means %>% select(-y), 
            by = c("node", "group")) %>% 
  group_by(node, group) %>% 
  summarise(final_mu = mean(final_mu), 
            y = mean(y))

# rpart
m0 <- rpart::rpart(y ~ X1, data_k1$train)
m1 <- rpart::rpart(y ~ X1 + group, data_k1$train)

rmse_results[9, 4] <- rss(data_k1$test$y, predict(m0, data_k1$test))
rmse_results[9, 5] <- rss(data_k1$test$y, predict(m1, data_k1$test))

# rf
m0 <- ranger::ranger(y ~ X1, data_k1$train, num.trees = 100)
m1 <- ranger::ranger(y ~ X1 + group, data_k1$train, num.trees = 100)
rmse_results[9, 6] <- rss(data_k1$test$y, predict(m0,  data_k1$test)$predictions)
rmse_results[9, 7] <- rss(data_k1$test$y, predict(m1,  data_k1$test)$predictions)
rmse_results$model <- as.character(rmse_results$model)
rmse_results$model[9] <- "2 node, sampling k1 (bigger)"
#-----------------------------------------------------------------------

rmse_results <- rmse_results %>% 
  mutate_at(2:7, ~round(as.numeric(.x), 3))
rmse_results$model.data <- rmse_results$model
rmse_results <- rmse_results %>% select(-model)
names(rmse_results)[1] <- "model/data"

write.table(rmse_results, file = "results/rmse_results.txt")

print(rmse_results %>% 
        xtable::xtable(), include.rownames = FALSE)
