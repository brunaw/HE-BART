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
rmse <- function(x, y) sqrt(mean((x - y)**2))
#-------------------------------------------------------------------------------
# Predictions are now looking weird; bug in the code?
set.seed(20)
data_sim <- sim_bart_add(n = 500, m = 10, k2 = 5, k1 = 8)
data_sim$tau # tau = 1.67

group_variable = "group"
formula <- y ~ X1
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 5, alpha = alpha, beta = beta, mu_mu = 0
)

data_sim$data$y <- c(scale(data_sim$data$y))

data_sim$data %>% 
  ggplot(aes(y = y, x = factor(group), group = group)) +
  geom_boxplot(fill = "#c95a49", alpha = 0.7) +
  theme_light(16) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(x = 'Group') 
  ggsave("paper/sim_data.png")

data_split <- data_sim$data %>% 
  mutate(split = ifelse(runif(n()) < 0.75, "train", "test")) %>% 
  split(.$split)

train <- data_split$train %>% select(-split)
test <- data_split$test %>% select(-split)
# -------------------------------------------------------------------------
m0_2node <- bart(formula, 
                 dataset = train,
                 iter = 100, 
                 P = 20, 
                 group_variable, pars, 
                 scale = FALSE)

pred_hbart <- predict_hbm(m0_2node, test, 
                          formula, group_variable)
pred_hbart_train <- predict_hbm(m0_2node, train, 
                          formula, group_variable)

# LME3  ------------------------------------------------------------------------
# 85% of data -------
lm3_m0_normal  <- lmer(y ~ X1 + (1 |group), data = train)
pred_lm3       <- predict(lm3_m0_normal, test)
pred_lm3_train <- predict(lm3_m0_normal, train)
rmse(pred_lm3, test$y)
# Bayesian LME  ----------------------------------------------------------------
pr = prior(normal(0, 1), class = 'b')

# 85% of data -------
#blme <-  brm(
#   y ~ X1 + (1 |group), data = train, prior = pr, cores = 4)
# 
# pred_blme       <- predict(blme, test)
# pred_blme_train <- predict(blme, train)
pred_blme <- 0
pred_blme_train <- 0

# BART  ------------------------------------------------------------------------
# 85% of data -------
bart_0 <-  dbarts::bart2(y ~ X1, 
                         data = train,
                         test = test, 
                         keepTrees = TRUE)
pred_bart_train <- bart_0$yhat.train.mean
pred_bart <- bart_0$yhat.test.mean


all_preds <- data.frame(
  y = c(test$y, train$y),
  pred_hbart = c(pred_hbart$pred, pred_hbart_train$pred), 
  pred_lme = c(pred_lm3, pred_lm3_train),
  #pred_blme = c(pred_blme[, 1], pred_blme_train[, 1]), 
  pred_blme =  0, 
  pred_bart = c(pred_bart, pred_bart_train), 
  source = rep(c("Test (25%)", "Train (75%)"), c(nrow(test), nrow(train)))
) %>% 
  mutate(
    res_hbart = y - pred_hbart,
    res_lme = y - pred_lme,
    res_blme = y - pred_blme,
    res_bart = y - pred_bart
  ) %>% 
  dplyr::select(6:10) 
  #group_by(source) %>% 
  #summarise_all(~round(sqrt(mean(.x^2)), 3))
  #mutate(scenario = "All IDs in training data") %>% 

all_preds %>% 
group_by(source) %>%
summarise_all(~round(sqrt(mean(.x^2)), 3))

# source      res_hbart res_lme res_blme res_bart
# Test (25%)      0.916   0.941    1.04      1.01
# Train (75%)     0.7     0.88     0.984     0.89


all_wrangle <- all_preds %>% 
  pivot_longer(cols = c(res_hbart,res_lme, res_blme, res_bart))  %>% 
  mutate(name = as.factor(name))

levels(all_wrangle$name) <- c("BART", "BLME", 
                              "HE-BART", "LME")

all_wrangle %>% 
  ggplot(aes(y = value, x = name)) +
  geom_boxplot(fill = "#c95a49", alpha = 0.7) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  facet_wrap(~source, ncol = 2) +
  labs(y = "Residuals", x = "Algorithm") +
  theme_light(16) 
ggsave("paper/sim_residuals.png", 
       width = 7, height = 4)

#------------------------------------------------------------
df_tau <- data.frame(tau = m0_2node$tau_post) %>%
  slice(-1) %>%
  mutate(ind = 1:n()) %>% 
  slice(-c(1:15))

true_tau <- tau
p22 <- ggplot(df_tau, aes(
  x = 1/sqrt(tau)
)) +
  geom_vline(xintercept = 1/sqrt(true_tau),
             colour = '#c95a49', size = 1.2, linetype = 'dotted') +
  geom_density(alpha = 0.4) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x = expression('Sampled values of 1/sqrt('~tau~')'), y = "Density",
       title = expression('Dotted line = true 1/sqrt('~tau~")")) +
  theme_bw(10)
p22
ggsave("paper/sim_tau.png", 
       width = 3, height = 2)

df_k1 <- data.frame(k1 = m0_2node$sampled_k1) %>%
  slice(-1) %>%
  mutate(ind = 1:n()) %>% 
  slice(-c(1:15))

true_k1 <- k1
p11 <- ggplot(df_k1, aes(
  x = k1
)) +
  geom_vline(xintercept = true_k1,
             colour = '#c95a49', size = 1.2, linetype = 'dotted') +
  geom_density(alpha = 0.4) +
  labs(x = expression('Sampled values of '~k[1]), y = "Density",
       title = expression('Dotted line = true '~k[1])) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_bw(10)
p11
ggsave("paper/sim_k1.png", 
       width = 3, height = 2)
#------------------------------------------------------------

