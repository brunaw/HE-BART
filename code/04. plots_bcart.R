#-----------------------------------------------------------------
all_data <- readRDS("data/all_data.rds")
data_one_node <- readRDS("data/data_one_node.rds")
data_k1 <- readRDS("data/data_k1.rds")
library(tidyverse)
files <- list.files("code/mixedbart/R") %>% 
  paste0("code/mixedbart/R/", .)
map(c(files[-c(7, 10, 11)]), source)
iter <-  3500
lim <- 500
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
# > tau
# [1] 0.8339095
tau <- 1.29

res_2node <- data.frame(
  taus = m0_2node$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m0_2node$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m0_2node$mus, ~{.x[2]})[-c(1:lim)]
)

res_2node %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/tau_2node.png")

mu = c(1.377732, -1.438701)


res_2node %>% 
  select(2:3) %>% 
  gather() %>% 
  mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key) + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  geom_vline(aes(xintercept = true_mean), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/mu_2node.png")


m0_2node$final_tree %>%
  group_by(node) %>% 
  summarise(m = max(X1), 
            mean_mu = mean(mu_sampled),
            mean_mu_j = mean(mu_js_sampled)) 

final_tree <- m0_2node$final_tree %>%
  group_by(node) %>% 
  mutate(m = as.factor(round(max(X1), 1)))

true_mus <- expand_grid(group = c(1:5),
                        m = unique(final_tree$m)) %>%
  arrange(m, group) %>% 
  mutate(muj = c(muj_2, muj_1))

final_tree <- final_tree %>% 
  right_join(true_mus, by = c("group", "m"))

final_tree %>% 
  ggplot(aes(x = y)) +
  geom_density() +
  geom_vline(aes(xintercept = mu_js_sampled), 
             colour = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = muj), colour = "blue") +
  facet_wrap(m~group, nrow = 2, scales = "free_x") +
  theme_bw()
ggsave(file = "results/ys_2node.png", height = 5, width = 10)
#-----------------------------------------------------------------
#-----------------------------------------------------------------
tau <- 2
# 2. 3 node tree, small ks
# res_3node <- data.frame(
#   taus = m0_3node$tau_post[-c(1:lim)],
#   mu_1 = map_dbl(m0_3node$mus, ~{.x[1]})[-c(1:lim)],
#   mu_2 = map_dbl(m0_3node$mus, ~{.x[2]})[-c(1:lim)]
# )
# 
# res_3node %>% 
#   ggplot(aes(x = sqrt(1/taus))) +
#   geom_density(fill = "grey90") +
#   geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
#              colour = "red", linetype = "dashed", size = 1) +
#   geom_vline(aes(xintercept = sqrt(1/tau)), 
#              colour = "blue", size = 1) +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) + 
#   theme_bw()
# ggsave(file = "results/tau_3node.png")
# # -------------------------------------------------------------------
# mujs <- c(sim_ks_small$mujs[[1]], 
#           sim_ks_small$mujs[[3]], 
#           sim_ks_small$mujs[[2]])
# 
# tree_mus <- m0_3node$final_tree %>%
#   group_by(node, group) %>% 
#   summarise(
#     mean_mu_j = mean(mu_js_sampled), 
#     mean_y = mean(y)) %>% 
#   ungroup() %>% 
#   mutate(mean_real = mujs) %>% 
#   left_join(m0_3node$final_tree, by = c("node", "group"))
# 
# tree_mus %>% 
#   group_by(node, group) %>% 
#   ggplot(aes(x = y)) +
#   geom_density() +
#   facet_wrap(node~group, nrow = 3, scales = "free_x") +
#   geom_vline(aes(xintercept = mean_mu_j), 
#              colour = "red", linetype = "dashed") +
#   geom_vline(aes(xintercept = mean_real), colour = "blue") +
#   theme_bw()
# 
# ggsave(file = "results/ys_mus_2node_tree.png", 
#        height = 10, width = 10)
# 
# 
# res_bcart %>% 
#   select(2:4) %>% 
#   gather() %>% 
#   mutate(true_mean = rep(sim_ks_small$mu[c(1,2,3)], 
#                          each = (iter-lim+1))) %>%
#   group_by(key) %>% 
#   mutate(mean_sample = mean(value)) %>% 
#   ggplot(aes(x = value)) +
#   geom_density(fill = "grey90") +
#   facet_wrap(~key) + 
#   geom_vline(aes(xintercept = mean_sample), 
#              colour = "red", linetype = "dashed") + 
#   geom_vline(aes(xintercept = true_mean), 
#              colour = "blue", size = 1) +
#   theme_bw()
# ggsave(file = "results/mus_2node_tree.png", height = 5, width = 10)
# -------------------------------------------------------------------
# 2. 3 node tree, small ks
res_3node <- data.frame(
  taus = m0_3node$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m0_3node$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m0_3node$mus, ~{.x[2]})[-c(1:lim)],
  mu_3 = map_dbl(m0_3node$mus, ~{.x[3]})[-c(1:lim)]
)

res_3node %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) + 
  theme_bw()
ggsave(file = "results/tau_3node.png")
# -------------------------------------------------------------------
mujs <- c(sim_ks_small$mujs[[1]], 
          sim_ks_small$mujs[[2]], 
          sim_ks_small$mujs[[3]])

tree_mus <- m0_3node$final_tree %>%
  group_by(node, group) %>% 
  summarise(
    mean_mu_j = mean(mu_js_sampled), 
    mean_y = mean(y)) %>% 
  ungroup() %>% 
  mutate(mean_real = mujs) %>% 
  left_join(m0_3node$final_tree, by = c("node", "group"))

tree_mus %>% 
  group_by(node, group) %>% 
  ggplot(aes(x = y)) +
  geom_density() +
  facet_wrap(node~group, nrow = 3, scales = "free_x") +
  geom_vline(aes(xintercept = mean_mu_j), 
             colour = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = mean_real), colour = "blue") +
  theme_bw()

ggsave(file = "results/ys_mujs_3node.png", 
       height = 10, width = 10)

res_3node %>% 
  select(2:4) %>% 
  gather() %>% 
  mutate(true_mean = rep(sim_ks_small$mu[c(1,2,3)], 
                         each = (iter-lim+1))) %>%
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key) + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  geom_vline(aes(xintercept = true_mean), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/mus_3node.png", height = 5, width = 10)
# -------------------------------------------------------------------
# 2. 3 node tree, big k2 
res_3node2 <- data.frame(
  taus = m1_3node$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m1_3node$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m1_3node$mus, ~{.x[2]})[-c(1:lim)],
  mu_3 = map_dbl(m1_3node$mus, ~{.x[3]})[-c(1:lim)]
)

res_3node2 %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) + 
  theme_bw()
ggsave(file = "results/tau_3node_big_k2.png")
# -------------------------------------------------------------------
mujs <- c(sim_k2_big$mujs[[1]], 
          sim_k2_big$mujs[[2]], 
          sim_k2_big$mujs[[3]])

tree_mus <- m1_3node$final_tree %>%
  group_by(node, group) %>% 
  summarise(
    mean_mu_j = mean(mu_js_sampled), 
    mean_y = mean(y)) %>% 
  ungroup() %>% 
  mutate(mean_real = mujs) %>% 
  left_join(m1_3node$final_tree, by = c("node", "group"))

tree_mus %>% 
  group_by(node, group) %>% 
  ggplot(aes(x = y)) +
  geom_density() +
  facet_wrap(node~group, nrow = 3, scales = "free_x") +
  geom_vline(aes(xintercept = mean_mu_j), 
             colour = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = mean_real), colour = "blue") +
  theme_bw()

ggsave(file = "results/ys_mus_3node_big_k2.png", 
       height = 10, width = 10)

res_3node2 %>% 
  select(2:4) %>% 
  gather() %>% 
  mutate(true_mean = rep(sim_k2_big$mu[c(1,2,3)], 
                         each = (iter-lim+1))) %>%
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key) + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  geom_vline(aes(xintercept = true_mean), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/mus_3node_big_k2.png", height = 5, width = 10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# 3. 3 node tree, big k1
res_3node3 <- data.frame(
  taus = m2_3node$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m2_3node$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m2_3node$mus, ~{.x[2]})[-c(1:lim)],
  mu_3 = map_dbl(m2_3node$mus, ~{.x[3]})[-c(1:lim)]
)

res_3node3 %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) + 
  theme_bw()
ggsave(file = "results/tau_3node_big_k1.png")
# -------------------------------------------------------------------
mujs <- c(sim_k1_big$mujs[[1]], 
          sim_k1_big$mujs[[2]], 
          sim_k1_big$mujs[[3]])

tree_mus <- m2_3node$final_tree %>%
  group_by(node, group) %>% 
  summarise(
    mean_mu_j = mean(mu_js_sampled), 
    mean_y = mean(y)) %>% 
  ungroup() %>% 
  mutate(mean_real = mujs) %>% 
  left_join(m2_3node$final_tree, by = c("node", "group"))

tree_mus %>% 
  group_by(node, group) %>% 
  ggplot(aes(x = y)) +
  geom_density() +
  facet_wrap(node~group, nrow = 3, scales = "free_x") +
  geom_vline(aes(xintercept = mean_mu_j), 
             colour = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = mean_real), colour = "blue") +
  theme_bw()

ggsave(file = "results/ys_mus_3node_big_k1.png", 
       height = 10, width = 10)

res_3node3 %>% 
  select(2:4) %>% 
  gather() %>% 
  mutate(true_mean = rep(sim_k1_big$mu[c(2,1,3)], 
                         each = (iter-lim+1))) %>%
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key) + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  geom_vline(aes(xintercept = true_mean), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/mus_3node_big_k1.png", height = 5, width = 10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# 4. 3 node tree, big ks
res_3node4 <- data.frame(
  taus = m3_3node$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m3_3node$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m3_3node$mus, ~{.x[2]})[-c(1:lim)],
  mu_3 = map_dbl(m3_3node$mus, ~{.x[3]})[-c(1:lim)]
)

res_3node4 %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) + 
  theme_bw()
ggsave(file = "results/tau_3node_big_ks.png")
# -------------------------------------------------------------------
mujs <- c(sim_ks_big$mujs[[1]], 
          sim_ks_big$mujs[[2]], 
          sim_ks_big$mujs[[3]])

tree_mus <- m3_3node$final_tree %>%
  group_by(node, group) %>% 
  summarise(
    mean_mu_j = mean(mu_js_sampled), 
    mean_y = mean(y)) %>% 
  ungroup() %>% 
  mutate(mean_real = mujs) %>% 
  left_join(m3_3node$final_tree, by = c("node", "group"))

tree_mus %>% 
  group_by(node, group) %>% 
  ggplot(aes(x = y)) +
  geom_density() +
  facet_wrap(node~group, nrow = 3, scales = "free_x") +
  geom_vline(aes(xintercept = mean_mu_j), 
             colour = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = mean_real), colour = "blue") +
  theme_bw()

ggsave(file = "results/ys_mus_3node_big_ks.png", 
       height = 10, width = 10)

res_3node4 %>% 
  select(2:4) %>% 
  gather() %>% 
  mutate(true_mean = rep(sim_ks_big$mu[c(1,2,3)], 
                         each = (iter-lim+1))) %>%
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key) + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  geom_vline(aes(xintercept = true_mean), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/mus_3node_big_ks.png", height = 5, width = 10)
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# 5. Fitting a regular tree

m0_regular$final_tree$node %>% n_distinct()
res_regular <- data.frame(
  taus = m0_regular$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m0_regular$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m0_regular$mus, ~{.x[2]})[-c(1:lim)]
)

res_regular %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) + 
  theme_bw()
ggsave(file = "results/tau_regular.png")
# -------------------------------------------------------------------
tree_mus <- m0_regular$final_tree %>%
  group_by(node, group) %>% 
  mutate(
    mean_mu_j = mean(mu_js_sampled), 
    mean_y = mean(y)) 

m0_regular$final_tree %>%
  group_by(node, group) %>% 
  summarise(
    mean_mu_j = mean(mu_js_sampled), 
    mean_y = mean(y)) 

tree_mus %>% 
  group_by(node, group) %>% 
  ggplot(aes(x = y)) +
  geom_density() +
  geom_rug() + 
  facet_wrap(node~group, nrow = 2, scales = "free_x") +
  geom_vline(aes(xintercept = mean_mu_j), 
             colour = "red", linetype = "dashed") +
  theme_bw()

ggsave(file = "results/ys_mus_regular.png", 
       height = 7, width = 10)

res_regular %>% 
  select(2:3) %>% 
  gather() %>% 
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key, scales = "free_x") + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  theme_bw()
ggsave(file = "results/mus_regular.png", height = 4, width = 6)
# -------------------------------------------------------------------
#----------------------------------------------------------------- 
# 6. Sampling k1, fixed one-split tree
# > tau
# [1] 1.535821
tau <- 1.535821
iter <- 750
lim <- 150

res_sample_comp <- data.frame(
  taus = m0_k1_comp$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m0_k1_comp$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m0_k1_comp$mus, ~{.x[2]})[-c(1:lim)]
)
res_sample_comp %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  theme_bw()

ggsave(file = "results/tau_samp_k1_comp.png")

mu = c(-1.628563, 1.483883)

res_sample_comp %>% 
  select(2:3) %>% 
  gather() %>% 
  mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key) + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  geom_vline(aes(xintercept = true_mean), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/mu_samp_k1_comp.png")


m0_k1_comp$final_tree %>%
  group_by(node) %>% 
  summarise(m = max(X1), 
            mean_mu = mean(mu_sampled),
            mean_mu_j = mean(mu_js_sampled)) 

final_tree <- m0_k1_comp$final_tree %>%
  group_by(node) %>% 
  mutate(m = as.factor(round(max(X1), 1)))

true_mus <- expand_grid(group = c(1:25),
                        m = unique(final_tree$m)) %>%
  arrange(m, group) %>% 
  mutate(muj = c(muj_2, muj_1))


final_tree <- final_tree %>% 
  right_join(true_mus, by = c("group", "m"))

final_tree %>% 
  ggplot(aes(x = y)) +
  geom_density() +
  geom_vline(aes(xintercept = mu_js_sampled), 
             colour = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = muj), colour = "blue") +
  facet_wrap(m~group, nrow = 5, scales = "free_x") +
  theme_bw()

ggsave(file = "results/ys_samp_k1_comp.png",
       height = 12, width = 10)

# mean((m0_k1_comp$final_tree$y - 
#         m0_k1_comp$final_tree$mu_js_sampled)^2)
# mean((m0_k1_samp$final_tree$y - 
#         m0_k1_samp$final_tree$mu_js_sampled)^2)

res_sample_k1 <- data.frame(
  taus = m0_k1_samp$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m0_k1_samp$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m0_k1_samp$mus, ~{.x[2]})[-c(1:lim)],
  k1 = m0_k1_samp$sampled_k1[-c(2:lim)]
)

res_sample_k1 %>% 
  group_by(k1) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(mean_k1 = mean(k1)) %>% 
  ggplot(aes(x = k1)) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean_k1), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 8), 
             colour = "blue", size = 1) +
  theme_bw()

ggsave(file = "results/k1_samples.png")


res_sample_k1 %>% 
  # group_by(k1) %>%
  # slice(1) %>%
  ungroup() %>% 
  mutate(mean_k1 = mean(k1), iteration = 1:n()) %>% 
  ggplot(aes(y = k1, x = iteration)) +
  geom_point(fill = "grey90") +
  geom_hline(aes(yintercept = mean_k1), 
             colour = "red", linetype = "dashed", size = 1) +
  theme_bw()

ggsave(file = "results/k1_samples_chain.png")


res_sample_k1 %>% 
  group_by(k1) %>%
  slice(1) %>%
  ungroup() %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  theme_bw()

ggsave(file = "results/tau_samp_k1.png")

res_sample_k1 %>% 
  group_by(k1) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(mean_k1_tau = mean(k1/taus)) %>% 
  ggplot(aes(x = sqrt(k1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = sqrt(mean_k1_tau)), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(k1_tau)), 
             colour = "blue", size = 1) +
  #xlim(0, 2) +
  theme_bw()
ggsave(file = "results/k1_tau.png")

mu = c(-1.628563, 1.483883)

res_sample_k1 %>% 
  select(2:3) %>% 
  gather() %>% 
  mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key) + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  geom_vline(aes(xintercept = true_mean), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/mu_samp_k1.png")


m0_k1_samp$final_tree %>%
  group_by(node) %>% 
  summarise(m = max(X1), 
            mean_mu = mean(mu_sampled),
            mean_mu_j = mean(mu_js_sampled)) 

final_tree <- m0_k1_samp$final_tree %>%
  group_by(node) %>% 
  mutate(m = as.factor(round(max(X1), 1)))

true_mus <- expand_grid(group = c(1:25),
                        m = unique(final_tree$m)) %>%
  arrange(m, group) %>% 
  mutate(muj = c(muj_2, muj_1))


final_tree <- final_tree %>% 
  right_join(true_mus, by = c("group", "m"))

final_tree %>% 
  ggplot(aes(x = y)) +
  geom_density() +
  geom_vline(aes(xintercept = mu_js_sampled), 
             colour = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = muj), colour = "blue") +
  facet_wrap(m~group, nrow = 5, scales = "free_x") +
  theme_bw()
ggsave(file = "results/ys_samp_k1.png", height = 12, width = 10)



# res_sample_k1 %>% 
#   ggplot(aes(x = k1)) +
#   geom_density() +
#   geom_vline(aes(xintercept = 7), 
#              colour = "red", linetype = "dashed") +
#   theme_bw()
# ggsave(file = "results/ys_2node.png", height = 5, width = 10)
#-----------------------------------------------------------------
# 7. Sampling k1, fixed one-split tree, bigger unif
# > tau
# [1] 1.535821
tau <- 1.535821
iter <- 750
lim <- 150

k1_tau <- 8/tau


res_sample_k1_big <- data.frame(
  taus = m0_k1_samp_bigger$tau_post[-c(1:lim)],
  mu_1 = map_dbl(m0_k1_samp_bigger$mus, ~{.x[1]})[-c(1:lim)],
  mu_2 = map_dbl(m0_k1_samp_bigger$mus, ~{.x[2]})[-c(1:lim)],
  k1 = m0_k1_samp_bigger$sampled_k1[-c(2:lim)]
)

res_sample_k1_big %>% 
  group_by(k1) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(mean_k1 = mean(k1)) %>% 
  ggplot(aes(x = k1)) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean_k1), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = 8), 
             colour = "blue", size = 1) +
  theme_bw()

ggsave(file = "results/k1_samples_big.png")

res_sample_k1_big %>% 
  # group_by(k1) %>%
  # slice(1) %>%
  ungroup() %>% 
  mutate(mean_k1 = mean(k1), iteration = 1:n()) %>% 
  ggplot(aes(y = k1, x = iteration)) +
  geom_point(fill = "grey90") +
  geom_hline(aes(yintercept = mean_k1), 
             colour = "red", linetype = "dashed", size = 1) +
  theme_bw()
ggsave(file = "results/k1_samples_chain_big.png")


res_sample_k1_big %>% 
  ggplot(aes(x = sqrt(1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = mean(sqrt(1/taus))), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(1/tau)), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/tau_samp_k1_big.png")


res_sample_k1_big %>% 
  group_by(k1) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(mean_k1_tau = mean(k1/taus)) %>% 
  ggplot(aes(x = sqrt(k1/taus))) +
  geom_density(fill = "grey90") +
  geom_vline(aes(xintercept = sqrt(mean_k1_tau)), 
             colour = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = sqrt(k1_tau)), 
             colour = "blue", size = 1) +
  #xlim(0, 2) +
  theme_bw()
ggsave(file = "results/k1_tau_big.png")

c(muj_1, muj_2) %>% sd()


mu = c(-1.628563, 1.483883)

res_sample_k1_big %>% 
  select(2:3) %>% 
  gather() %>% 
  mutate(true_mean = rep(mu[c(2, 1)], each = (iter-lim+1))) %>%
  group_by(key) %>% 
  mutate(mean_sample = mean(value)) %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey90") +
  facet_wrap(~key) + 
  geom_vline(aes(xintercept = mean_sample), 
             colour = "red", linetype = "dashed") + 
  geom_vline(aes(xintercept = true_mean), 
             colour = "blue", size = 1) +
  theme_bw()
ggsave(file = "results/mu_samp_k1_big.png")


m0_k1_samp_bigger$final_tree %>%
  group_by(node) %>% 
  summarise(m = max(X1), 
            mean_mu = mean(mu_sampled),
            mean_mu_j = mean(mu_js_sampled)) 

final_tree <- m0_k1_samp_bigger$final_tree %>%
  group_by(node) %>% 
  mutate(m = as.factor(round(max(X1), 1)))

true_mus <- expand_grid(group = c(1:25),
                        m = unique(final_tree$m)) %>%
  arrange(m, group) %>% 
  mutate(muj = c(muj_2, muj_1))


final_tree <- final_tree %>% 
  right_join(true_mus, by = c("group", "m"))

final_tree %>% 
  ggplot(aes(x = y)) +
  geom_density() +
  geom_vline(aes(xintercept = mu_js_sampled), 
             colour = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = muj), colour = "blue") +
  facet_wrap(m~group, nrow = 5, scales = "free_x") +
  theme_bw()

ggsave(file = "results/ys_samp_k1_big.png", 
       height = 12, width = 10)


# res_sample_k1 %>% 
#   ggplot(aes(x = k1)) +
#   geom_density() +
#   geom_vline(aes(xintercept = 7), 
#              colour = "red", linetype = "dashed") +
#   theme_bw()
# ggsave(file = "results/ys_2node.png", height = 5, width = 10)
#-----------------------------------------------------------------