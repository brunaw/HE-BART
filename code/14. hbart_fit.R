#-------------------------------------------------------------------------------
# A Hierarchical BART model run for a one-split, fixed tree  
# 
#-------------------------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(ranger)
library(lme4)
library(brms)
library(dbarts)

source("code/00. simulation_functions.R")

files <- list.files("code/hmb/R") %>% 
  paste0("code/hmb/R/", .)
map(files, source)

rss <- function(x, y) sum((x - y)**2)
#-------------------------------------------------------------------------------
# Predictions are now looking weird; bug in the code?
set.seed(15)
data_sim <- sim_bart(n = 1000, m = 10, k2 = 5, k1 = 8)
data_sim$tau # tau = 1.67

group_variable = "group"
formula <- y ~ X1
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 5, alpha = alpha, beta = beta, mu_mu = 0
)

data_sim$data$y <- c(scale(data_sim$data$y))

data_split <- data_sim$data %>% 
  mutate(split = ifelse(runif(n()) < 0.75, "train", "test")) %>% 
  split(.$split)

train <- data_split$train %>% select(-split)
test <- data_split$test %>% select(-split)
# -------------------------------------------------------------------------
m0_2node <- bart(formula, 
                 dataset = train,
                 iter = 100, 
                 P = 50, 
                 group_variable, pars, 
                 scale = FALSE)

mean(m0_2node$sampled_k1)
mean(m0_2node$tau_post)

# My BART version ---------------------------------------------------------
pars$k1 <- 0.001

m0_2node_bart <- bart(formula, 
                      dataset = train,
                      iter = 75, 
                      P = 15,
                      group_variable, pars, 
                      scale = FALSE, sample_k1 = FALSE)


# -------------------------------------------------------------------------
pred_m0 <- predict_hbm(m0_2node, newdata = test, formula = formula, 
                       group_variable = group_variable)

rss(test$y, pred_m0$pred)
pred_m0_bart <- predict_hbm(m0_2node_bart, newdata = test, formula = formula,
                            group_variable = group_variable)
rss(test$y, pred_m0_bart$pred)

# rr <- pred_m0 %>%
#   pivot_longer(cols = c(y, pred)) %>%
#   mutate(name = as.factor(name))
# 
# levels(rr$name) <- c("Prediction", "True y")
# rr %>%
#   ggplot(aes(value, y = factor(group), group = group)) +
#   geom_boxplot(colour = 'green4') +
#   facet_wrap(~ name) +
#   labs(x = 'Values', y = 'Groups') +
#   guides(colour = 'none') +
#   theme_bw()


# cc_bart <- m0_2node_bart$final_trees %>% 
#   mutate(current_tree = map(tree_data, ~{ tail(.x, 1)})) %>% 
#   select(tree_index, current_tree) %>% 
#   unnest(current_tree) %>% 
#   select(tree_index, est_tree) %>% 
#   unnest(est_tree) %>% 
#   dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>% 
#   ungroup()


# res_bart <- cc_bart %>% 
#   dplyr::select(tree_index, y, sampled_mu_j, group) %>% 
#   pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>% 
#   mutate(pred = rowSums(.[3:ncol(.)]), 
#          res = (y - pred)**2, 
#          res_2 = (y - pred)) 

# -------------------------------------------------------------------------
cc <- m0_2node$final_trees %>%
  mutate(current_tree = map(tree_data, ~{ tail(.x, 1)})) %>%
  select(tree_index, current_tree) %>%
  unnest(current_tree) %>%
  select(tree_index, est_tree) %>%
  unnest(est_tree) %>%
  dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>%
  ungroup()


res <- cc %>%
  dplyr::select(tree_index, y, sampled_mu_j, group) %>%
  pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>%
  mutate(pred = rowSums(.[3:ncol(.)]),
         res = (y - pred)**2,
         res_2 = (y - pred))
# 

# -------------------------------------------------------------------------
# Tau traceplot
df_tau <- data.frame(tau = m0_2node$tau_post) %>%
  slice(-1) %>%
  mutate(ind = 1:n())

tau_TRUE <- data_sim$tau

res_sd <- sd(res$res_2)
1/sqrt(mean(m0_2node$tau_post))
sq_tau_true <- 1/sqrt(tau_TRUE)


p22 <- ggplot(df_tau, aes(
  x = 1/sqrt(tau)
)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = sq_tau_true, colour = 'Green4') +
  geom_vline(xintercept = res_sd, colour = 'Red') +
  labs(x = expression('Sampled values of 1/sqrt('~tau~')'), y = "Density",
       title = 'Green line = true 1/sqrt(tau);\n Red line = residual SD') +
  theme_bw()
p22

p23 <-  ggplot(df_tau, aes(
  x = ind, y = 1/sqrt(tau), 
)) +
  geom_point() +
  geom_hline(aes(yintercept = median(1/sqrt(tau))), colour = "red") +
  labs(y = expression('Sampled values of 1/sqrt('~tau~')'), x = "Iteration") +
  theme_bw()


p22 + p23 + plot_layout(nrow = 2)
ggsave(file = "results/reports/february/tau.png")

df_k1 <- data.frame(k1 = m0_2node$sampled_k1) %>%
  slice(-1) %>%
  group_by(k1) %>%
  slice(1)
pars$k1 <- 8

p33 <- ggplot(df_k1, aes(
  x = k1
)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = pars$k1, colour = 'Green4') +
  labs(x = expression('Sampled values of '~k[1]), y = "Density") +
  theme_bw()
p33


p44 <- data.frame(k1 = m0_2node$sampled_k1)  %>% 
  bind_cols(df_tau) %>% 
  ggplot(aes(
  y = k1/tau, x = ind
)) +
  geom_point() +
  geom_hline(aes(yintercept = median(k1/tau)), colour = "red") +
  labs(y = expression("k"[1]~'/'~tau), x = "Iteration") +
  theme_bw()
p44

# ggsave(file = "results/reports/february/k1_over_tau.png")
# 
#----------------------------------------------------------------------
# rr <- res %>% 
#   select(pred, y, group) %>% 
#   pivot_longer(cols = c(y, pred)) %>% 
#   mutate(name = as.factor(name)) 
# 
# levels(rr$name) <- c("Prediction", "True y")
# rr %>% 
#   ggplot(aes(value, y = factor(group), group = group)) +
#   geom_boxplot(colour = 'green4') +
#   facet_wrap(~ name) + 
#   labs(x = 'Values', y = 'Groups') +
#   guides(colour = 'none') +
#   theme_bw()
# 
# #ggsave(file = "results/reports/january/predictions.png")
# 
# res_mu <- cc %>%
#   dplyr::select(tree_index, y, sampled_mu, group) %>%
#   pivot_wider(names_from = tree_index, values_from = sampled_mu) %>%
#   mutate(pred = rowSums(.[3:ncol(.)]),
#          res = (y - pred)**2)
# # 
# # 
# # res_mu %>% 
# #   ggplot(aes(y, pred)) +
# #   geom_point(aes(colour = factor(group)))
# # 
# # 
# # rr <- res_mu %>% 
# #   select(pred, y, group) %>% 
# #   pivot_longer(cols = c(y, pred)) %>% 
# #   mutate(name = as.factor(name))
# 
# 
# # rr %>% 
# #   ggplot(aes(value, y = factor(group), group = group)) +
# #   geom_boxplot(colour = 'green4') +
# #   facet_wrap(~ name) + 
# #   labs(x = 'Values', y = 'Groups') +
# #   guides(colour = 'none') +
# #   theme_bw()
# #ggsave(file = "presentations/group_meeting_2021/predictions.png")
# 
# #ggsave(file = "results/reports/december/predictions_fixed_one_node.png")

#-------------------------------------------------------------------------
lm3_m0 <- lmer(y ~ X1 + (1 |group), data = train)
pred_lm3_m0 <- predict(lm3_m0, test)
rss(pred_lm3_m0, test$y)

# Bayesian lm
pr = prior(normal(0, 1), class = 'b')

bayesian_mixed = brm(
  y ~ X1 + (1 |group) , 
  data = train,
  prior = pr,
  cores = 4
)

pred_bayesian_mixed <- predict(bayesian_mixed, test)

# BART -------------------------------------------------------------------------
bart_0 = dbarts::bart2(formula, 
                       data = train,
                       test = test, 
                       keepTrees = TRUE)
pred_bart <- bart_0$yhat.test.mean

# Also do the RSS per group
all_preds <- data.frame(
  group = test$group, 
  y = test$y, 
  pred_hbm = pred_m0$pred, 
  pred_my_bart = pred_m0_bart$pred, 
  pred_lme = c(pred_lm3_m0), 
  pred_blme = pred_bayesian_mixed[, 1], 
  pred_bart = pred_bart
) %>% 
  mutate(
    res_hbm = y - pred_hbm,
    res_my_bart = y - pred_my_bart, 
    res_lme = y - pred_lme,
    res_blme = y - pred_blme,
    res_bart = y - pred_bart,
      
  )

ssqrs <- all_preds %>% 
  select(8:12) %>% 
  summarise_all(~round(sum(.x^2))) %>% 
  pivot_longer(cols = c(res_hbm, res_my_bart, res_lme, res_blme, res_bart)) %>% 
  set_names(c("name", "ssqrs")) %>% 
  mutate(ssqrs_text = paste0("SSR: ", ssqrs))


all_wrangle <- all_preds %>% 
  pivot_longer(cols = c(res_hbm, res_my_bart,res_lme,  res_blme, res_bart)) %>% 
  mutate(name = as.factor(name)) %>% 
  left_join(ssqrs, by = 'name') %>% 
  select(name, value, ssqrs, ssqrs_text) %>% 
  mutate(name = as.factor(name))


levels(all_wrangle$name) <- c("BART", "Bayesian LME", 
                              "Hierachical BART", "LME", 
                              "My BART\n (k1 = 0)")

all_wrangle %>% 
  group_by(name) %>% 
  mutate(max = max(value)) %>% 
  ggplot(aes(y = value,  x = name, group = factor(name))) +
  geom_boxplot(colour = 'green4', aes(fill = ssqrs)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#32CD32") + 
  geom_text(aes(x=name, y=max+0.5, label = ssqrs_text), 
            col='black', size=7, position = position_dodge(width = .75)) +
  labs(x = 'Models', y = 'Residuals', title = "Sum of squared residuals in the test set") +
  guides(colour = 'none', fill = 'none') +
  theme_bw(16)

#ggsave(file = "results/reports/february/residuals_test.png")
# OK
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
df_real <- sleepstudy %>% 
  set_names(c('y', 'X1', 'group'))

df_real$y <- c(scale(df_real$y))

df_real %>% 
  ggplot(aes(y = y, x = factor(group),  group = factor(group))) +
  geom_boxplot() +
  theme_bw(18) +
  labs(x = 'Group')

#ggsave(file = "results/reports/february/sleepstudy.png")

data_split <- df_real %>% 
  group_by(group) %>% 
  mutate(split = ifelse(runif(n()) < 0.85, "train", "test")) %>% 
  split(.$split)

train <- data_split$train %>% select(-split) %>% ungroup()
test <- data_split$test %>% select(-split) %>% ungroup()
table(train$group)
table(test$group)
formula <- y ~ X1 

# HBM --------------------------------------------------------------
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 5, alpha = alpha, beta = beta, mu_mu = 0
)
group_variable <- "group"
m0_real <- bart(formula, 
                dataset = train,
                iter = 150, 
                P = 50, 
                group_variable, pars, 
                min_u = 0, max_u = 20, 
                scale = FALSE)

# My BART --------------------------------------------------------------
pars$k1 <- 0.001

m0_real_bart <- bart(formula, 
                dataset = train,
                iter = 150, 
                P = 50, 
                group_variable, pars, 
                scale = FALSE, sample_k1 = FALSE)

# -------------------------------------------------------------------------
pred_m0 <- predict_hbm(model = m0_real, newdata = test, formula = formula,
                       group_variable = group_variable)

pred_m0_bart <- predict_hbm(model = m0_real_bart, newdata = test, formula = formula,
                            group_variable = group_variable)
# # -----------------------------------------------------------------------
# cc <- m0_real$final_trees %>%
#   mutate(current_tree = map(tree_data, ~{ tail(.x, 1)})) %>%
#   select(tree_index, current_tree) %>%
#   unnest(current_tree) %>%
#   select(tree_index, est_tree) %>%
#   unnest(est_tree) %>%
#   dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>%
#   ungroup()
# 
# res <- cc %>%
#   dplyr::select(tree_index, y, sampled_mu_j, group) %>%
#   pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>%
#   mutate(pred = rowSums(.[3:ncol(.)]),
#          res = (y - pred)**2,
#          res_2 = (y - pred))
# 
# cc <- m0_real_bart$final_trees %>%
#   mutate(current_tree = map(tree_data, ~{ tail(.x, 1)})) %>%
#   select(tree_index, current_tree) %>%
#   unnest(current_tree) %>%
#   select(tree_index, est_tree) %>%
#   unnest(est_tree) %>%
#   dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>%
#   ungroup()
# 
# res_bart <- cc %>%
#   dplyr::select(tree_index, y, sampled_mu_j, group) %>%
#   pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>%
#   mutate(pred = rowSums(.[3:ncol(.)]),
#          res = (y - pred)**2,
#          res_2 = (y - pred))


lm3_m0 <- lmer(y ~ X1 + (1 |group), data = train)
pred_lm3_m0 <- predict(lm3_m0, test)
pred_lm3_m0_train <- predict(lm3_m0)


# Bayesian lm
pr = prior(normal(0, 1), class = 'b')

bayesian_mixed = brm(
  y ~  X1 + (1 |group) , 
  data = train,
  prior = pr,
  cores = 4
)

pred_bayesian_mixed <- predict(bayesian_mixed, test)
pred_bayesian_mixed_train <- predict(bayesian_mixed)


# Bart
bart_0 = dbarts::bart2(formula, data = train, test = test)
pred_bart <- bart_0$yhat.test.mean
pred_bart_train <- bart_0$yhat.train.mean


# TEST ssr
all_preds <- data.frame(
  group = test$group, 
  y = test$y, 
  pred_hbm = pred_m0$pred, 
  pred_my_bart = pred_m0_bart$pred, 
  pred_lme = c(pred_lm3_m0), 
  pred_blme = pred_bayesian_mixed[, 1], 
  pred_bart = pred_bart
) %>% 
  mutate(
    res_hbm = y - pred_hbm,
    res_my_bart = y - pred_my_bart, 
    res_lme = y - pred_lme,
    res_blme = y - pred_blme,
    res_bart = y - pred_bart,
    
  )

ssqrs <- all_preds %>% 
  select(8:12) %>% 
  summarise_all(~round(sum(.x^2), 1)) %>% 
  pivot_longer(cols = c(res_hbm, res_my_bart, res_lme, res_blme, res_bart)) %>% 
  set_names(c("name", "ssqrs")) %>% 
  mutate(ssqrs_text = paste0("SSR: ", ssqrs))


all_wrangle <- all_preds %>% 
  pivot_longer(cols = c(res_hbm, res_my_bart,res_lme,  res_blme, res_bart)) %>% 
  mutate(name = as.factor(name)) %>% 
  left_join(ssqrs, by = 'name') %>% 
  select(name, value, ssqrs, ssqrs_text) %>% 
  mutate(name = as.factor(name))


levels(all_wrangle$name) <- c("BART", "Bayesian LME", 
                              "Hierachical BART", "LME", 
                              "My BART\n (k1 = 0)")

all_wrangle %>% 
  group_by(name) %>% 
  mutate(max = max(value)) %>% 
  ggplot(aes(y = value,  x = name, group = factor(name))) +
  geom_boxplot(colour = 'green4', aes(fill = ssqrs)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#32CD32") + 
  geom_text(aes(x=name, y=max+0.5, label = ssqrs_text), 
            col='black', size=7, position = position_dodge(width = .75)) +
  labs(x = 'Models', y = 'Residuals', title = "Sum of squared residuals in the test set") +
  guides(colour = 'none', fill = 'none') +
  theme_bw(18)

#ggsave(file = "results/reports/february/residuals_sleepstudy_test.png")

# Train -----------------------------------------------------------------------
# Train ssr
# all_preds <- data.frame(
#   group = train$group, 
#   y = train$y, 
#   pred_hbm = res$pred, 
#   pred_my_bart = res_bart$pred, 
#   pred_lme = c(pred_lm3_m0_train), 
#   pred_blme = pred_bayesian_mixed_train[, 1], 
#   pred_bart = pred_bart_train
# ) %>% 
#   mutate(
#     res_hbm = y - pred_hbm,
#     res_my_bart = y - pred_my_bart, 
#     res_lme = y - pred_lme,
#     res_blme = y - pred_blme,
#     res_bart = y - pred_bart,
#     
#   )
# 
# ssqrs <- all_preds %>% 
#   select(8:12) %>% 
#   summarise_all(~round(sum(.x^2), 1)) %>% 
#   pivot_longer(cols = c(res_hbm, res_my_bart, res_lme, res_blme, res_bart)) %>% 
#   set_names(c("name", "ssqrs")) %>% 
#   mutate(ssqrs_text = paste0("SSR: ", ssqrs))
# 
# 
# all_wrangle <- all_preds %>% 
#   pivot_longer(cols = c(res_hbm, res_my_bart,res_lme,  res_blme, res_bart)) %>% 
#   mutate(name = as.factor(name)) %>% 
#   left_join(ssqrs, by = 'name') %>% 
#   select(name, value, ssqrs, ssqrs_text) %>% 
#   mutate(name = as.factor(name))
# 
# 
# levels(all_wrangle$name) <- c("BART", "Bayesian LME", 
#                               "Hierachical BART", "LME", 
#                               "My BART\n (k1 = 0)")
# 
# all_wrangle %>% 
#   group_by(name) %>% 
#   mutate(max = max(value)) %>% 
#   ggplot(aes(y = value,  x = name, group = factor(name))) +
#   geom_boxplot(colour = 'green4', aes(fill = ssqrs)) +
#   scale_fill_gradient(low = "#FFFFFF", high = "#32CD32") + 
#   geom_text(aes(x=name, y=max+0.5, label = ssqrs_text), 
#             col='black', size=7, position = position_dodge(width = .75)) +
#   labs(x = 'Models', y = 'Residuals', title = "Sum of squared residuals in the train set") +
#   guides(colour = 'none', fill = 'none') +
#   theme_bw(18)
# 
# ggsave(file = "results/reports/ferbuary/residuals_sleepstudy_train.png")
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# wages is not a very good dataset for trees

# library(brolgar)
# data(wages)
# 
# alpha = 0.5; beta = 1; mu_mu = 0;
# pars <- list(
#   k1 = 8, k2 = 5, alpha = alpha, beta = beta, mu_mu = 0
# )
# 
# 
# first_150 <- unique(wages$id)[1:150]
# wages <- wages %>% filter(id %in% first_150)
# 
# filter_out <- wages %>% 
#   count(id) %>% 
#   arrange(n) %>% 
#   filter(n < 4)
# 
# wages <- wages %>% filter(!id %in% filter_out$id)
# dim(wages)
# 
# names(wages)[2] <- "y"
# table(wages$xp)
# formula <- y ~  black + unemploy_rate + ged +  high_grade  + xp_since_ged + hispanic
# group_variable <-  "id"
# 
# test <-   wages %>%
#   group_by(id) %>% 
#   slice(1)
# 
# train <- wages %>% 
#   group_by(id) %>% 
#   slice(2:n())
# 
# nrow(train)/nrow(wages)
# nrow(test)/nrow(wages)
# 
# # # Model -------------------------------------------------------------
# m0_real <- bart(formula,
#                 dataset = train,
#                 iter = 100,
#                 P = 50,
#                 group_variable, pars,
#                 min_u = 0, max_u = 20,
#                 sample_k1 = TRUE)
#  
#    
# # # My BART --------------------------------------------------------------
# pars$k1 <- 0.001
# 
# m0_real_bart <- bart(formula,
#                      dataset = train,
#                      iter = 25,
#                      P = 15,
#                      group_variable, pars,
#                      scale = FALSE, sample_k1 = FALSE)
# 
# # -------------------------------------------------------------------------
# cc <- m0_real$final_trees %>%
#   mutate(current_tree = map(tree_data, ~{ tail(.x, 1)})) %>%
#   select(tree_index, current_tree) %>%
#   unnest(current_tree) %>%
#   select(tree_index, est_tree) %>%
#   unnest(est_tree) %>%
#   dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>%
#   ungroup()
# 
# 
# res <- cc %>%
#   dplyr::select(tree_index, y, sampled_mu_j, group) %>%
#   pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>%
#   mutate(pred = rowSums(.[3:ncol(.)]),
#          res = (y - pred)**2,
#          res_2 = (y - pred))
# 
# # Tau traceplot
# df_tau <- data.frame(tau = m0_real$tau_post) %>%
#   slice(-1) %>%
#   mutate(ind = 1:n())
# 
# #tau_TRUE <- data_sim$tau
# res_sd <- sd(res$res_2)
# #1/sqrt(mean(m0_2node$tau_post))
# #sq_tau_true <- 1/sqrt(tau_TRUE)
# 
# p22 <- ggplot(df_tau, aes(
#   x = 1/sqrt(tau)
# )) +
#   geom_vline(xintercept = res_sd, colour = 'Red') +
#   geom_density(alpha = 0.4) +
#   labs(x = expression('Sampled values of 1/sqrt('~tau~')'), y = "Density",
#        title = 'Red line = residual SD') +
#   theme_bw()
# p22
# 
# p23 <-  ggplot(df_tau, aes(
#   x = ind, y = 1/sqrt(tau), 
# )) +
#   geom_point() +
#   geom_hline(aes(yintercept = median(1/sqrt(tau))), colour = "red") +
#   labs(y = expression('Sampled values of 1/sqrt('~tau~')'), x = "Iteration") +
#   theme_bw()
# 
# 
# p22 + p23 + plot_layout(nrow = 2)
# # ggsave(file = "results/reports/february/tau_wages.png")
# 
# 
# # #-------------------------------------------------------------------------
# pred_m0 <- predict_hbm(model = m0_real, newdata = test, formula = formula,
#                        group_variable = group_variable)
# 
# pred_m0_bart <- predict_hbm(
#   model = m0_real_bart, newdata = test, formula = formula,
#   group_variable = group_variable)
# # 
# # #-------------------------------------------------------------------------
# lm3_m0 <- lmer(y ~ black + unemploy_rate + ged +  high_grade + xp + (1 |id), data = train)
# pred_lm3_m0 <- predict(lm3_m0, test)
# rss(test$y, pred_m0$pred)
# rss(test$y, pred_lm3_m0)
# 
# bayesian_mixed = brm(
#   y ~ black + unemploy_rate + ged +  high_grade + xp + (1 |id) ,
#   data = train,
#   prior = pr,
#   cores = 4
# )
#  
# pred_bayesian_mixed <- predict(bayesian_mixed, test)
# 
# # BART -------------------------------------------------------------------------
# bart_0 = dbarts::bart2(formula,
#                        data = train,
#                        test = test,
#                        keepTrees = TRUE)
# pred_bart <- bart_0$yhat.test.mean
# 
# # -------------------------------------------------------------------------
# # Common splits?
# all_results <- select(m0_real$final_trees, results) %>% 
#   mutate(id = 1:n()) %>% 
#   unnest(results) %>% 
#   drop_na()
# 
# all_results2 <- select(m0_real_bart$final_trees, results) %>% 
#   mutate(id = 1:n()) %>% 
#   unnest(results) %>% 
#   drop_na()
# 
# cart <- rpart::rpart(formula, train)
# pp <- predict(cart, test)
# rss(pp, test$y)
# # No splits at all ------ how does this compare to BART?
# # n.trees = 75L, n.samples = 500L, n.burn = 500L, nchains = 4
# all_results
# bart_trees <- bart_0$fit$getTrees(75)
# wages %>% summary()
# bart_trees %>% 
#   View()
# max(bart_trees$sample)
# table(bart_trees$chain)
# bart_0$fit$printTrees(1)
# 
# # -------------------------------------------------------------------------
# all_preds <- data.frame(
#   group = test$id,
#   y = test$y,
#   pred_hbm = pred_m0$pred,
#   pred_my_bart = pred_m0_bart$pred,
#   pred_lme = c(pred_lm3_m0),
#   pred_blme = pred_bayesian_mixed[, 1],
#   pred_bart = pred_bart
# ) %>%
#   mutate(
#     res_hbm = y - pred_hbm,
#     res_my_bart = y - pred_my_bart,
#     res_lme = y - pred_lme,
#     res_blme = y - pred_blme,
#     res_bart = y - pred_bart,
# 
#   )
# 
# ssqrs <- all_preds %>%
#   select(8:12) %>%
#   summarise_all(~round(sum(.x^2), 2)) %>%
#   pivot_longer(cols = c(res_hbm, res_my_bart, res_lme, res_blme, res_bart)) %>%
#   set_names(c("name", "ssqrs")) %>%
#   mutate(ssqrs_text = paste0("SSR: ", ssqrs))
# 
# 
# all_wrangle <- all_preds %>%
#   pivot_longer(cols = c(res_hbm, res_my_bart,res_lme,  res_blme, res_bart)) %>%
#   mutate(name = as.factor(name)) %>%
#   left_join(ssqrs, by = 'name') %>%
#   select(name, value, ssqrs, ssqrs_text) %>%
#   mutate(name = as.factor(name))
# 
# 
# levels(all_wrangle$name) <- c("BART", "Bayesian LME",
#                               "Hierachical BART", "LME",
#                               "My BART\n (k1 = 0)")
# 
# all_wrangle %>%
#   group_by(name) %>%
#   mutate(max = max(value)) %>%
#   ggplot(aes(y = value,  x = name, group = factor(name))) +
#   geom_boxplot(colour = 'green4', aes(fill = ssqrs)) +
#   scale_fill_gradient(low = "#FFFFFF", high = "#32CD32") +
#   geom_text(aes(x=name, y=max+0.5, label = ssqrs_text),
#             col='black', size=7, position = position_dodge(width = .75)) +
#   labs(x = 'Models', y = 'Residuals', title = "Sum of squared residuals in the test set") +
#   guides(colour = 'none', fill = 'none') +
#   theme_bw(16)
# 
# ggsave(file = "results/reports/february/residuals_test_wages.png")
# # 
# # 
# # 
# # 
# # #-------------------------------------------------------------------------
# # # Friedman data
# # df <- sim_friedman_bart(n = 500, j = 10)
# # 
# # pars <- list(
# #   k1 = 8, k2 = 10, alpha = alpha, beta = beta, mu_mu = 0
# # )
# # 
# # # Bug: predictions are not varying
# # df %>% 
# #   ggplot(aes(y = y, x = factor(group),  group = factor(group))) +
# #   geom_boxplot() +
# #   theme_bw(18) +
# #   labs(x = 'Group')
# # 
# # #ggsave(file = "results/reports/january/sim_fried.png")
# # 
# # formula <- y ~ X1 + X2 + X3 + X4 + X5
# # df$y <- c(scale(df$y))
# # m0_fried <- bart(formula, 
# #                  dataset = df,
# #                  iter = 50, 
# #                  P = 50, 
# #                  group_variable, pars, 
# #                  scale = FALSE)
# # m0_fried$sampled_k1 %>% mean()
# # 
# # cc <- m0_fried$final_trees %>% 
# #   mutate(current_tree = map(tree_data, ~{ slice(.x, (50))})) %>% 
# #   select(tree_index, current_tree) %>% 
# #   unnest(current_tree) %>% 
# #   select(tree_index, est_tree) %>% 
# #   unnest(est_tree) %>% 
# #   dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>% 
# #   ungroup()
# # 
# # 
# # res <- cc %>% 
# #   dplyr::select(tree_index, y, sampled_mu_j, group) %>% 
# #   pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>% 
# #   mutate(pred = rowSums(.[3:ncol(.)]), 
# #          res = (y - pred)**2, 
# #          res_2 = (y - pred)) 
# # 
# # 
# # rr <- res %>% 
# #   select(pred, y, group) %>% 
# #   pivot_longer(cols = c(y, pred)) %>% 
# #   mutate(name = as.factor(name)) 
# # 
# # 
# # levels(rr$name) <- c("Prediction", "True y")
# # 
# # 
# # rr %>% 
# #   ggplot(aes(value, y = factor(group), group = group)) +
# #   geom_boxplot(colour = 'green4') +
# #   facet_wrap(~ name) + 
# #   labs(x = 'Values', y = 'Groups') +
# #   guides(colour = 'none') +
# #   theme_bw()
# # 
# # #ggsave(file = "results/reports/january/predictions_fried.png")
# # 
# # lm3_m0 <- lmer(y ~ X1 + X2 + X3  + X4 + X5 + (1 |group), data = df)
# # pred_lm3_m0 <- predict(lm3_m0)
# # 
# # # Bayesin lm
# # pr = prior(normal(0, 1), class = 'b')
# # 
# # bayesian_mixed = brm(
# #   y ~  X1 + X2 + X3 +  X4 + X5 + (1 |group) , 
# #   data = df,
# #   prior = pr,
# #   cores = 4
# # )
# # 
# # pred_bayesian_mixed <- predict(bayesian_mixed)
# # 
# # x_mat <- df %>% 
# #   select(X1, X2, X3, X4, X5) %>% 
# #   as.matrix()
# # 
# # # Bart
# # bart_0 = dbarts::bart2(x_mat, df$y)
# # pred_bart <- bart_0$yhat.train.mean
# # 
# # 
# # all_preds <- data.frame(
# #   group = df$group, 
# #   y = df$y, 
# #   pred_hbm = res$pred, 
# #   pred_lme = c(pred_lm3_m0), 
# #   pred_blme = pred_bayesian_mixed[, 1], 
# #   pred_bart = pred_bart
# # ) %>% 
# #   mutate(
# #     res_hbm = y - pred_hbm,
# #     res_lme = y - pred_lme,
# #     res_blme = y - pred_blme,
# #     res_bart = y - pred_bart,
# #     
# #   )
# # 
# # cor(all_preds$y, all_preds$pred_hbm)
# # cor(all_preds$y, all_preds$pred_blme)
# # 
# # 
# # ssqrs <- all_preds %>% 
# #   select(7:10) %>% 
# #   summarise_all(~round(sum(.x^2, 2))) %>% 
# #   pivot_longer(cols = c(res_hbm, res_lme, res_blme, res_bart)) %>% 
# #   set_names(c("name", "ssqrs")) %>% 
# #   mutate(ssqrs = paste0("SSR = ", ssqrs))
# # 
# # 
# # all_wrangle <- all_preds %>% 
# #   pivot_longer(cols = c(res_hbm, res_lme, res_blme, res_bart)) %>% 
# #   mutate(name = as.factor(name)) %>% 
# #   left_join(ssqrs, by = 'name') %>% 
# #   select(name, value, ssqrs) %>% 
# #   mutate(name = as.factor(name))
# # 
# # levels(all_wrangle$name) <- c("BART", "Bayesian LME", 
# #                               "Hierachical BART", "LME")
# # all_wrangle %>% 
# #   group_by(name) %>% 
# #   mutate(max = max(value)) %>% 
# #   ggplot(aes(y = value,  x = name, group = factor(name))) +
# #   geom_boxplot(colour = 'green4') +
# #   geom_text(aes(x=name, y=max+1, label = ssqrs), 
# #             col='black', size=5, position = position_dodge(width = .75)) +
# #   labs(x = 'Models', y = 'Residuals') +
# #   guides(colour = 'none') +
# #   theme_bw(18)
# # #ggsave(file = "results/reports/january/residuals_sim_fried.png")
# # 
# # #-------------------------------------------------------------------------
# # 
# # 
# # #-------------------------------------------------------------------------
# # # --> https://www.rensvandeschoot.com/tutorials/lme4/
# # library(haven)
# # 
# # popular2data <- read_sav(
# #   file ="https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%202/popularity/SPSS/popular2.sav?raw=true")
# # 
# # 
# # popular2data <- select(
# #   popular2data, pupil, class, extrav, sex, texp, popular)
# # 
# # popular2data$sex <- as.factor(popular2data$sex)
# # n_distinct(popular2data$class) # 100 classes, could this work?
# # glimpse(popular2data)
# # 
# # 
# # popular2data %>% 
# #   ggplot(aes(y = popular, x = factor(class),  group = factor(class))) +
# #   geom_boxplot() +
# #   theme_bw(18) +
# #   labs(x = 'Group')
# # 
# # group_variable = "class"
# # formula <- popular ~ sex + extrav + texp
# # alpha = 0.5; beta = 1; mu_mu = 0;
# # pars <- list(
# #   k1 = 8, k2 = 5, alpha = alpha, beta = beta, mu_mu = 0
# # )
# # 
# # popular2data$popular <- c(scale(popular2data$popular))
# # m0_pop <- bart(formula, 
# #                dataset = popular2data, 
# #                iter = 10, 
# #                P = 50, 
# #                group_variable, pars, 
# #                scale = FALSE)
# # 
# # 
# # 
# # lm3_pop <- lmer(popular ~ 1 + sex + extrav + texp + (1 | class), 
# #                data = popular2data)
# # pred_lm3 <- predict(lm3_pop)
# # summary(model2)
# # 
# # 
# # #-------------------------------------------------------------------------------
# # # NOT IN USE
# # #-------------------------------------------------------------------------
# # 
# # 
# # bart_1 = dbarts::bart2(cbind(data_sim$data$X1, data_sim$data$group), data_sim$data$y)
# # cor(data_sim$data$y, bart_0$yhat.train.mean) 
# # cor(data_sim$data$y, bart_1$yhat.train.mean)
# # cor(res$pred, data_sim$data$y)
# # 
# # 
# # 
# # sum((c(bart_1$yhat.train.mean) - data_sim$data$y)^2)
# # sum(res$res)
# # sum(res_mu$res)
# # 
# # library(BART)
# # 
# # m0_bart <- gbart(data_sim$data$X1, 
# #               data_sim$data$y, ndpost = 500, ntree = 5)
# # m1_bart <- gbart(cbind(data_sim$data$X1, data_sim$data$group), 
# #                  data_sim$data$y, ndpost = 500, ntree = 5)
# # pp_0 <- m0_bart$yhat.train.mean
# # pp_1 <- m1_bart$yhat.train.mean
# # sum((c(pp_0) - data_sim$data$y)^2)
# # sum((c(pp_1) - data_sim$data$y)^2)
# # #pp <- predict(post, newdata = matrix(data_sim$data$X1))
# # 
# # # rss_df$rss[6] <- sum(res$res)
# # # rss_df$rss[5] <- sum(res_mu$res)
# # 
# # 
# # 
# # 
# # #------------------------------------------------------------------
# # # rpart
# # m0 <- rpart::rpart(y ~ X1, data_sim$data)
# # m1 <- rpart::rpart(y ~ X1 + group, data_sim$data)
# # rss_df$rss[1] <- rss(data_sim$data$y, predict(m0, data_sim$data))
# # rss_df$rss[2]  <- rss(data_sim$data$y, predict(m1, data_sim$data))
# # 
# # # rf
# # m0 <- ranger::ranger(y ~ X1, data_sim$data, num.trees = 5)
# # m1 <- ranger::ranger(y ~ X1 + group, data_sim$data, num.trees = 5)
# # rss_df$rss[3] <- rss(data_sim$data$y, predict(m0,  data_sim$data)$predictions)
# # rss_df$rss[4] <- rss(data_sim$data$y, predict(m1,  data_sim$data)$predictions)
# # rss_df$rss[6] <- sum(res$res)
# # 
# # 
# # rss_df$rss[5] <- sum(res_mu$res)
# # 
# # rss_df$source <- c("tree", "tree", "RF", "RF", "HBart", "HBart")
# # rss_df$type <- rep(c("normal", "group"), 3)
# # rss_df
# # # rss_df %>% 
# # #   select(source, type, rss) %>% 
# # #   write.table(file = "results/reports/december/rss.txt")
# # 
# # # > m0 <- rpart::rpart(y ~ X1, data_sim$data)
# # # > m1 <- rpart::rpart(y ~ X1 + group, data_sim$data)
# # # > rss(data_sim$data$y, predict(m0, data_sim$data))
# # # [1] 9201.931
# # # > rss(data_sim$data$y, predict(m1, data_sim$data))
# # # [1] 1209.953
# # # > # rf
# # #   > m0 <- ranger::ranger(y ~ X1, data_sim$data, num.trees = 5)
# # # > m1 <- ranger::ranger(y ~ X1 + group, data_sim$data, num.trees = 5)
# # # > rss(data_sim$data$y, predict(m0,  data_sim$data)$predictions)
# # # [1] 3925.085
# # # > rss(data_sim$data$y, predict(m1,  data_sim$data)$predictions)
# # # [1] 1225.62
# # # > sum(res$res)
# # # [1] 1175.315
# # 
# # # plot group densities (for all trees) and compare it to 
# # # true averages 
# # 
# # library(BART)
# # 
# # post <- gbart(cbind(data_sim$data$X1, data_sim$data$group), data_sim$data$y, ndpost = 1000)
# # pp <- post$yhat.train.mean
# # #pp <- predict(post, newdata = matrix(data_sim$data$X1))
# # sum((c(pp) - data_sim$data$y)^2)
# # 
# # 
# # # Calculate mean of all_trees
# # 
# # # calc_tree <- function(i, data = m0_2node$final_trees){
# # #   cc <- data %>% 
# # #     mutate(current_tree = map(tree_data, ~{ slice(.x, (i))})) %>% 
# # #     select(tree_index, current_tree) %>% 
# # #     unnest(current_tree) %>% 
# # #     select(tree_index, est_tree) %>% 
# # #     unnest(est_tree) %>% 
# # #     dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>% 
# # #     ungroup()
# # #   
# # #   
# # #   pred <- cc %>% 
# # #     dplyr::select(tree_index, y, sampled_mu_j, group) %>% 
# # #     pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>% 
# # #     mutate(pred = rowSums(.[3:ncol(.)]))  %>% 
# # #     select(y, group, pred) %>% 
# # #     mutate(n_tree = i, ind = 1:n())
# # #  pred 
# # # }
# # 
# # 
# # # all_res <- map_dfr(.x = 5:150, calc_tree)
# # # summarise_res <- all_res %>% 
# # #   group_by(ind) %>% 
# # #   summarise(y = mean(y), group = mean(group), 
# # #             pred = mean(pred))
# # # 
# # # sum((all_res$y - all_res$pred)^2)
# # # summary(all_res$pred)
# # # summary(all_res$y)
