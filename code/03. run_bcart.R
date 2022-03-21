#-----------------------------------------------------------------
all_data <- readRDS("data/all_data.rds")
data_k1 <- readRDS("data/data_k1.rds")
library(tidyverse)
files <- list.files("code/hmb/R") %>% 
  paste0("code/hmb/R/", .)

map(files[-5], source)

#-----------------------------------------------------------------
attach(all_data)
group_variable = "group"
formula <- y ~ X1
iter <-  3500
lim <- 500
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)

# 2-node fixed tree ------- 
m0_2node <- bcart_fixed_one(formula, 
                            dataset = data_one_node$train,
                            iter = iter, 
                            group_variable, pars, 
                            scale = FALSE)
saveRDS(m0_2node, "results/m0_2node.rds")
# 3-node fixed tree, small ks ------- 
group_variable = "group"
formula <- y ~ X1 + X2
pars <- list(
  k1 = 1, k2 = 3, alpha = alpha, beta = beta, mu_mu = 0
)

m0_3node <- bcart_fixed(formula, 
                        dataset = sim_ks_small$data$train,
                        iter = iter, 
                        group_variable, pars, scale = FALSE)
saveRDS(m0_3node, "results/m0_3node.rds")
# 3-node fixed tree, bigger k2 ------- 
pars$k2 <- 10
pars$k1 <- 1.5
m1_3node <- bcart_fixed(formula, dataset = sim_k2_big$data$train, 
                        iter = iter, 
                        group_variable, pars, scale = FALSE)

saveRDS(m1_3node, "results/m1_3node.rds")
# 3-node fixed tree, bigger k1 ------- 
pars$k2 <- 15
pars$k1 <- 2.5
m2_3node <- bcart_fixed(formula, 
                        dataset = sim_k1_big$data$train, 
                        iter = iter, 
                        group_variable, pars, scale = FALSE)
saveRDS(m2_3node, "results/m2_3node.rds")
# 3-node fixed tree, bigger k1 and k2 ------- 
pars$k2 <- 10
pars$k1 <- 13
m3_3node <- bcart_fixed(formula, 
                        dataset = sim_ks_big$data$train, 
                        iter = iter, 
                        group_variable, pars, scale = FALSE)
saveRDS(m3_3node, "results/m3_3node.rds")
#-----------------------------------------------------------------
# Regular tree with two splits ---------------------------------
pars$k2 <- 10
pars$k1 <- 13
formula <- y ~ X1 + X2

m0_regular <- bcart(formula, 
                    dataset = sim_ks_big$data$train, 
                    iter = iter, 
                    group_variable, pars, scale = FALSE)
saveRDS(m0_regular, "results/m0_regular.rds")
#-----------------------------------------------------------------
# 1 node fixed tree with sampling of k1
iter <- 750
formula <- y ~ X1
group_variable <- "group"
pars <- list(
  k1 = 6, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)

# 2-node fixed tree with k1 sampling
m0_k1_comp <- bcart_fixed_one(
  formula, 
  dataset = data_k1$train,
  iter = iter, 
  group_variable, pars, 
  scale = FALSE)
saveRDS(m0_k1_comp, "results/m0_k1_comp.rds")

# 2-node fixed tree ------- 
m0_k1_samp <- bcart_fixed_k1(
  formula, 
  dataset = data_k1$train,
  iter = iter, 
  group_variable, pars, 
  scale = FALSE, min_u = 5, max_u = 12)
saveRDS(m0_k1_samp, "results/m0_k1_samp.rds")
#-----------------------------------------------------------------
# 1 node fixed tree with sampling of k1 with a wider range for U()
iter <- 750
formula <- y ~ X1
group_variable <- "group"
pars <- list(
  k1 = 6, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
)

# 2-node fixed tree ------- 
m0_k1_samp <- bcart_fixed_k1(
  formula, 
  dataset = data_k1$train,
  iter = iter, 
  group_variable, pars, 
  scale = FALSE, 
  min_u = 2, max_u = 18)

saveRDS(m0_k1_samp, "results/m0_k1_bigger.rds")
#-----------------------------------------------------------------
# 1 node fixed tree with sampling of k1 with a wider range for U()
# Not being used atm
# iter <- 2000
# formula <- y ~ X1
# group_variable <- "group"
# pars <- list(
#   k1 = 6, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
# )
# 
# # 2-node fixed tree ------- 
# m0_k1_samp_0 <- bcart_fixed_k1(
#   formula, 
#   dataset = data_k1$train,
#   iter = iter, 
#   group_variable, pars, 
#   scale = FALSE, 
#   min_u = 0, max_u = 18)
# 
# current_tree = m0_k1_samp_0$final_tree
# m0_k1_samp_0$sampled_k1 %>% mean()
# saveRDS(m0_k1_samp, "results/m0_k1_0.rds")
#-----------------------------------------------------------------
