library(tidyverse)
library(nlme)
library(tidymodels)
library(mlmRev)
library(lme4)
data("Exam")


select <- dplyr::select

files <- list.files("code/hmb/R") %>% 
  paste0("code/hmb/R/", .)
map(files, source)

data <- Exam %>% 
  select(school, standLRT, normexam) %>% 
  set_names(c('group', 'X1', 'y'))

sp_data <- data %>% 
  group_by(group) %>% 
  summarise(m = abs(cor(y, X1)))
gg <- sp_data %>% 
  arrange(m) %>% 
  slice(1:8)

data <- filter(data, group %in% gg$group)
dim(data)

data %>% 
  ggplot(aes(X1, y)) +
  geom_point() +
  facet_wrap(~group)

group_variable = "group"
formula <- y ~ X1
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 5, alpha = alpha, beta = beta, mu_mu = 0
)
rmse <- function(x, y) sqrt(mean((x - y)**2))

set.seed(2025)
data_split <- initial_split(data, prop = 0.90)
test <- testing(data_split)
train <- training(data_split)


sets <- list(train = train, test = test)
# This takes over a minute to start -----
m0_real <- bart(formula,
                dataset = sets$train,
                iter = 250,
                P = 50,
                group_variable, pars,
                min_u = 0, max_u = 20,
                scale = FALSE)

pred_hbart       <- predict_hbm(m0_real, sets$test, formula, group_variable)
pred_hbart_train <- predict_hbm(m0_real, sets$test, formula, group_variable)

rmse(pred_hbart$pred, sets$test$y)

# LME3  ------------------------------------------------------------------------
# 85% of data -------
lm3_m0_normal  <- lmer(y ~ X1 + (1 |group), data = train)
pred_lm3       <- predict(lm3_m0_normal, sets$test, 
                          re.form=NA)
pred_lm3_train <- predict(lm3_m0_normal, train)
rmse(pred_lm3, sets$test$y)
# Bayesian LME  ----------------------------------------------------------------
pr = prior(normal(0, 1), class = 'b')

# 85% of data -------
blme <-  brm(
  y ~ X1 + (1 |group), data = train, prior = pr, cores = 4)

pred_blme       <- predict(blme, test)
pred_blme_train <- predict(blme, train)


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
  pred_blme = c(pred_blme[, 1], pred_blme_train[, 1]), 
  pred_bart = c(pred_bart, pred_bart_train), 
  source = rep(c("Test (20%)", "Train (80%)"), c(nrow(test), nrow(train)))
) %>% 
  mutate(
    res_hbart = y - pred_hbart,
    res_lme = y - pred_lme,
    res_blme = y - pred_blme,
    res_bart = y - pred_bart
  ) %>% 
  dplyr::select(6:10) %>%
  group_by(source) %>% 
  summarise_all(~round(sqrt(mean(.x^2)), 3)) %>% 
  mutate(scenario = "All IDs in training data") %>% 
  select(6, 1, 2:5)


all_preds
