#-------------------------------------------------------------------------------
library(tidyverse)
library(tidymodels)
library(ranger)
library(lme4)
library(brms)
library(dbarts)
library(merTools)

select <- dplyr::select

files <- list.files("code/hmb/R") %>% 
  paste0("code/hmb/R/", .)
map(files, source)

rss <- function(x, y) sum((x - y)**2)
rmse <- function(x, y) sqrt(mean((x - y)**2))
set.seed(2022)
#-------------------------------------------------------------------------------
# General parameters 
group_variable = "group"
formula <- y ~ X1
alpha = 0.5; beta = 1; mu_mu = 0;
pars <- list(
  k1 = 8, k2 = 5, alpha = alpha, beta = beta, mu_mu = 0
)

#-------------------------------------------------------------------------------
# Dataset split 
df_real <- sleepstudy %>% 
  set_names(c('y', 'X1', 'group'))

df_real$y <- c(scale(df_real$y))

df_real %>% 
  ggplot(aes(y = y, x = factor(group),  group = factor(group))) +
  geom_boxplot() +
  theme_bw(18) +
  labs(x = 'Group')



set.seed(2025)
data_split <- initial_split(df_real, prop = 0.85)
test <- testing(data_split)
train <- training(data_split)

n_distinct(train$group)
n_distinct(test$group)

sets <- list(train = train, test = test)

# write_rds(sets, "paper/sets.rds")
#-------------------------------------------------------------------------------
sets <- read_rds("paper/sets.rds")
train <- sets$train
test <- sets$test

# HBART --------------------------------------------------------------
#-------------------------------------------------------------------------------
# Using 85% of training data 
#-------------------------------------------------------------------------------
# m0_real <- bart(formula, 
#                 dataset = train,
#                 iter = 300, 
#                 P = 50, 
#                 group_variable, pars, 
#                 min_u = 0, max_u = 20, 
#                 scale = FALSE)

plot_model <- m0_real
#write_rds(plot_model, file = "paper/plot_model.rds")
#-------------------------------------------------------------------------------
# Missing IDs in the training set
#-------------------------------------------------------------------------------
train_2 <- train %>% filter(!group %in% c(308, 309, 351))
test_2 <- test %>% bind_rows(train %>% filter(group %in% c(308, 309, 351)))
unique(test_2$group) %in% unique(train_2$group)

# m0_real <- bart(formula, 
#                 dataset = train_2,
#                 iter = 200, 
#                 P = 50, 
#                 group_variable, pars, 
#                 min_u = 0, max_u = 20, 
#                 scale = FALSE)
plot_model_missing <- m0_real
#write_rds(plot_model_missing, file = "paper/plot_model_missing.rds")
#------------------------------------------------------------------------------- 
# Reading models back ----------------------------------------------------------
plot_model <- read_rds("paper/plot_model.rds")
plot_model_missing <- read_rds("paper/plot_model_missing.rds")
#------------------------------------------------------------------------------- 
# Predictions 
pred_hbart       <- predict_hbm(plot_model, test, formula, group_variable)
pred_hbart_train <- predict_hbm(plot_model, train, formula, group_variable)

pred_hbart_missing       <- predict_hbm(plot_model_missing, test_2, formula, group_variable)
pred_hbart_train_missing <- predict_hbm(plot_model_missing, train_2, formula, group_variable)
# LME3  ------------------------------------------------------------------------
# 85% of data -------
lm3_m0_normal  <- lmer(y ~ X1 + (1 |group), data = train)
pred_lm3       <- predict(lm3_m0_normal, test)
pred_lm3_train <- predict(lm3_m0_normal, train)

# Missing IDs -------
lm3_m0_missing         <- lmer(y ~ X1 + (1 |group), data = train_2)
pred_lm3_missing       <- predict(lm3_m0_missing, test_2, re.form=NA)
pred_lm3_missing_train <- predict(lm3_m0_missing, train_2)

# Bayesian LME  ----------------------------------------------------------------
pr = prior(normal(0, 1), class = 'b')

# 85% of data -------
blme <-  brm(
  y ~ X1 + (1 |group), data = train, prior = pr, cores = 4)

pred_blme       <- predict(blme, test)
pred_blme_train <- predict(blme, train)


# Missing IDs -------
blme_missing <-  brm(
  y ~ X1 + (1 |group), data = train_2, prior = pr, cores = 4)

pred_blme_missing       <- predict(blme_missing, test_2, allow_new_levels = TRUE)
pred_blme_train_missing <- predict(blme_missing, train_2)

# BART  ------------------------------------------------------------------------
# 85% of data -------
bart_0 <-  dbarts::bart2(y ~ X1, 
                       data = train,
                       test = test, 
                       keepTrees = TRUE)
pred_bart_train <- bart_0$yhat.train.mean
pred_bart <- bart_0$yhat.test.mean


# Missing IDs -------
bart_missing = dbarts::bart2(y ~ X1, 
                             data = train_2,
                             test = test_2, 
                             keepTrees = TRUE)
pred_bart_missing_train <- bart_missing$yhat.train.mean
pred_bart_missing <- bart_missing$yhat.test.mean

#-------------------------------------------------------------------------------
# All predictions --------------------------------------------------------------

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
  summarise_all(~round(sqrt(mean(.x^2)), 2)) %>% 
  mutate(scenario = "All IDs in training data") %>% 
  select(6, 1, 2:5)

all_preds_missing <- data.frame(
  y = c(test_2$y, train_2$y),
  #group = c(test_2$group, train_2$group), 
  pred_hbart = c(pred_hbart_missing$pred, pred_hbart_train_missing$pred), 
  pred_lme = c(pred_lm3_missing, pred_lm3_missing_train),
  pred_blme = c(pred_blme_missing[, 1], pred_blme_train_missing[, 1]), 
  pred_bart = c(pred_bart_missing, pred_bart_missing_train), 
  source = rep(c("Test (25%)", "Train (75%)"), c(nrow(test_2), nrow(train_2)))
) %>% 
  mutate(
    res_hbart = y - pred_hbart,
    res_lme = y - pred_lme,
    res_blme = y - pred_blme,
    res_bart = y - pred_bart
  ) %>% 
  # dplyr::select(2, 6:10) %>%
  # group_by(source, group) %>% 
  dplyr::select(6:10) %>%
  group_by(source) %>% 
  summarise_all(~round(sqrt(mean(.x^2)), 2)) %>% 
  mutate(scenario = "Missing IDs in training data") %>% 
  select(6, 1, 2:5)

print(bind_rows(
  all_preds,
  all_preds_missing
) %>% 
  xtable::xtable(), include.rownames = FALSE)

# -----------------------------------------------------------------------------
# Plots -----------------------------------------------------------------------
ids <- unique(train$group)
range_X1 <- sort(unique(
  c(seq(min(train$X1), max(train$X1), length = 250), 
    sleepstudy$Days)))

new_test <- expand.grid(group = ids, X1 = range_X1)
pred_test <- predict_hbm(model = plot_model, newdata = new_test, 
                         formula = formula, 
                         group_variable = group_variable)
pred_test <-  
  bind_cols(pred_test, 
            predictInterval(lm3_m0_normal, new_test)) %>% 
  filter(group %in% c(308, 309, 351)) 

pred_test_final <- pred_test %>%
  mutate(low_ci = qnorm(0.025, mean = pred, sd = sqrt(var_mu)),
         upp_ci = qnorm(0.975, mean = pred, sd = sqrt(var_mu))) %>%
  left_join(df_real, by = c("group", "X1")) %>%
  mutate(group = paste0("ID: ", group),
         X1 = round(X1, 1),
         X_int = cut(X1, -0.0001:9))  %>% 
  group_by(group, X_int) %>%
  arrange(group, X1) %>%
  mutate(id = 1:n()) %>%
  mutate(
    # low_ci = ifelse(!id %in% c(1, 7, 14, 21),  NA, low_ci),
    # upp_ci = ifelse(!id %in% c(1, 7, 14, 21), NA, upp_ci),
    low_ci = ifelse(!id == 14,  NA, low_ci),
    upp_ci = ifelse(!id == 14, NA, upp_ci),
  )


rss_calc <- test %>% 
  filter(group %in% c(308, 309, 351)) 

pred_test_filter <- predict_hbm(model = plot_model, newdata = rss_calc, 
                                formula = formula, 
                                group_variable = group_variable)
pred_lme <- predict(lm3_m0_normal, rss_calc)
rss_hbart <- round(rmse(pred_test_filter$pred, rss_calc$y), 2)
rss_lme <- round(rmse(pred_lme, rss_calc$y), 2)

pred_test_final %>% 
  ggplot(aes(x = X1, y = pred)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill = "grey", alpha = 0.3) + 
  #geom_ribbon(aes(ymin=low_ci, ymax=upp_ci), fill = "#F96209", alpha = 0.3) +
  geom_line(aes(colour = "#F96209"), size = 0.5) +
  geom_line(colour = "#F96209", size = 0.5) +
  geom_point(aes(x = X1, y = y, colour = "#75E6DA"), size = 0.75) + 
  geom_point(aes(x = X1, y = y), colour = "#75E6DA", size = 2) + 
  geom_line(aes(x = X1, y = fit), colour = "grey") + 
  #geom_ribbon(aes(ymin = lwr_gam, ymax = fit_gam), fill = "lightblue") + 
  geom_errorbar(aes(ymin = low_ci, ymax = upp_ci),
                position = position_dodge(width = 0.2),
                width = 0.5, colour = '#F96209') +
  #geom_smooth(colour = "grey80",  alpha = 0.2) +
  facet_wrap(~group, ncol = 3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(y = "Average response time (ms)", 
       x = 'Covariate: days of sleep deprivation', 
       title = paste0("RMSE\nHE-BART: ", rss_hbart, ", LME: ", rss_lme)
  ) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_linedraw(20) +
  scale_colour_manual(
    name="Source:",
    values=c(Data="#75E6DA", 
             `HE-BART Prediction`="#F96209", 
              `LME Fit`= 'grey'), 
    guide = guide_legend(override.aes = list(
      size = c(3, 3, 3), shape = c(16, 16, 16)))) + 
  theme(panel.spacing.x = unit(0.5, "lines"), 
        legend.position = "bottom")

# ggsave(file = "results/reports/february/predictions_plot.png",
#        width = 12, height = 6)

ggsave(file = "paper/predictions_plot.png",
       width = 8, height = 5)

#-------------------------------------------------------------------------------
# -----------------------------------------------------------------------------

pred_test <- predict_hbm(model = plot_model_missing, newdata = new_test, 
                         formula = formula, 
                         group_variable = group_variable) %>% 
  filter(group %in% c(308, 309, 351)) 
new_test <- new_test %>% 
  filter(group %in% c(308, 309, 351)) 
# LME without a few groups is just lm? 
lm_missing <- lm(y ~ X1, data = train_2)

pred_test$fit <- predict(lm3_m0_missing, new_test, re.form=NA)
pred_lm <- predict(lm_missing, new_test, interval = "predict")
pred_test$upr <- pred_lm[, 3]
pred_test$lwr <- pred_lm[, 2]

pred_test_final <- pred_test %>%
  mutate(low_ci = qnorm(0.025, mean = pred, sd = sqrt(var_mu)),
         upp_ci = qnorm(0.975, mean = pred, sd = sqrt(var_mu))) %>%
  left_join(df_real, by = c("group", "X1")) %>%
  mutate(group = paste0("ID: ", group),
         X1 = round(X1, 1),
         X_int = cut(X1, -0.0001:9)) %>%
  group_by(group, X_int) %>%
  arrange(group, X1) %>%
  mutate(id = 1:n()) %>%
  mutate(
    low_ci = ifelse(id != 14, NA, low_ci),
    upp_ci = ifelse(id != 14, NA, upp_ci),
  )


rss_calc <- test_2 %>% 
  filter(group %in% c(308, 309, 351)) 

pred_test_filter <- predict_hbm(model = plot_model, newdata = rss_calc, 
                                formula = formula, 
                                group_variable = group_variable)
pred_lme <- predict(lm3_m0_normal, rss_calc)
rss_hbart <- round(rmse(pred_test_filter$pred, rss_calc$y), 2)
rss_lme <- round(rmse(pred_lme, rss_calc$y), 2)


pred_test_final %>% 
  ggplot(aes(x = X1, y = pred)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill = "grey", alpha = 0.3) + 
  geom_line(aes(colour = "#F96209"), size = 0.5) +
  geom_line(colour = "#F96209", size = 0.5) +
  geom_point(aes(x = X1, y = y, colour = "#75E6DA"), size = 0.75) + 
  geom_point(aes(x = X1, y = y), colour = "#75E6DA", size = 2) + 
  geom_line(aes(x = X1, y = fit), colour = "grey") + 
  #geom_ribbon(aes(ymin = lwr_gam, ymax = fit_gam), fill = "lightblue") + 
  geom_errorbar(aes(ymin = low_ci, ymax = upp_ci),
                position = position_dodge(width = 0.2),
                width = 0.5, colour = '#F96209') +
  #geom_smooth(colour = "grey80",  alpha = 0.2) +
  facet_wrap(~group, ncol = 3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(y = "Average response time (ms)", 
       x = 'Covariate: days of sleep deprivation', 
       title = paste0("RMSE\nHE-BART: ", rss_hbart, ", LME: ", rss_lme)
  ) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_linedraw(20) +
  scale_colour_manual(
    name="Source:",
    values=c(Data="#75E6DA", 
             `HE-BART Prediction`="#F96209", 
             `LME Fit`= 'grey'), 
    guide = guide_legend(override.aes = list(
      size = c(3, 3, 3), shape = c(16, 16, 16)))) + 
  theme(panel.spacing.x = unit(0.5, "lines"), 
        legend.position = "bottom")

ggsave(file = "paper/predictions_plot_missing.png",
       width = 8, height = 5)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------








#m0_real <- plot_model
# LM3 --------------------------------------------------------------------------


pred <- predict_hbm(plot_model, test, formula, group_variable)

# sum((pred$pred - test$y)^2)
# sum((pred_lm3_m0 - test$y)^2)


# Re-check changes in these prediction; should not happen [ALL OKAY NOW]
pred_lm3_m0_train <- predict(lm3_m0_normal, train)
pred_train <- predict_hbm(m0_real, train, formula, group_variable)


# data.frame(y = test$y, pred = pred$pred, lm3 = pred_lm3_m0) %>% 
#   rowwise() %>% 
#   mutate(rss = rss(y, pred), 
#          rsslm = rss(y, lm3)) %>% 
#   ungroup() %>% 
#   mutate(m = sum(rss), mm = sum(rsslm)) %>% 
#   print(n = 24)

rss_table2 <- rss_table %>% 
  add_row(
  lm3 = rss(pred_lm3_m0_train, train$y), 
  hbm = rss(pred_train$pred, train$y), 
  trained_in = "85% for training data", 
  type = "training RSS"
) %>% 
  add_row(
    lm3 = rss(pred_lm3_m0, test$y), 
    hbm = rss(pred$pred, test$y), 
    trained_in = "85% for training data", 
    type = "testing RSS"
  )

rss_table3 %>% 
  mutate_at(vars(1:2), ~round(.x, 1)) %>% 
  knitr::kable()

#write_rds(rss_table3, file = "paper/rss_table.rds")
# RSS sometimes depends on the sampled valued for the observation -- check that
# Using data that didn't have all IDs

pred <- predict_hbm(plot_model_missing, test_2, formula, group_variable)

# sum((pred_lm3_m0 - test_2$y)^2)
# sum((pred$pred - test_2$y)^2)

# Re-check changes in these prediction; should not happen [ALL OKAY NOW]
pred_lm3_m0_train <- predict(lm3_m0, train_2)
pred_train <- predict_hbm(m0_real, train_2, formula, group_variable)


rss_table3 <- rss_table2 %>% 
  add_row(
    lm3 = rss(pred_lm3_m0_train, train_2$y), 
    hbm = rss(pred_train$pred, train_2$y), 
    trained_in = "85% for training data, missing IDS", 
    type = "training RSS"
  ) %>% 
  add_row(
    lm3 = rss(pred_lm3_m0, test_2$y), 
    hbm = rss(pred$pred, test_2$y), 
    trained_in = "85% for training data, missing IDS", 
    type = "testing RSS"
  )

test_3 <- test_2 %>% filter(!group %in% c(308, 309, 351))

pred_lm3_m0 <- predict(lm3_m0, test_3, re.form=NA )
pred <- predict_hbm(m0_real, test_3, formula, group_variable)
rss(pred_lm3_m0, test_3$y)
rss(pred$pred, test_3$y)


# Using the latest model to plot
#-------------------------------------------------------------------------------
ids <- unique(train$group)
range_X1 <- sort(unique(
  c(seq(min(train$X1), max(train$X1), length = 250), 
    sleepstudy$Days)))

new_test <- expand.grid(group = ids, X1 = range_X1)
pred_test <- predict_hbm(model = plot_model, newdata = new_test, 
                         formula = formula, 
                         group_variable = group_variable)
pred_test <-  
  bind_cols(pred_test, 
  predictInterval(lm3_m0_normal, new_test)) %>% 
  filter(group %in% c(308, 309, 351)) 


# split_g <- pred_test %>% split(.$group)
# split_g <- list(split_g$`308`, split_g$`309`, split_g$`351`)
# 
# pred_split <- new_test %>% 
#   filter(group %in% c(308, 309, 351)) %>% 
#   split(.$group)
# pred_split <- list(pred_split$`308`, pred_split$`309`, pred_split$`351`)
# 
# fit_gam <- function(train_gam, pred_gam, group){
#   m <- gam(pred ~ s(X1), data = train_gam) 
#   p <- predict(m, newdata = pred_gam, se = TRUE) 
# 
#   g <- data.frame(X1 = pred_gam$X1, 
#                   fit_gam = p$fit,
#                   lwr_gam = p$fit - 5.96*p$se.fit, 
#                   upr_gam = p$fit + 5.96*p$se.fit, 
#                   group = group)
#   g
# }
# 
# groups <- c(308, 309, 351)
# gg <- pmap_dfr(list(split_g, pred_split, groups), fit_gam)
# gg$group <- as.factor(gg$group)

# pred_test2 <- pred_test %>% 
#   left_join(gg, by = c("X1", "group"))

pred_test_final <- pred_test %>%
  mutate(low_ci = qnorm(0.025, mean = pred, sd = sqrt(var_mu)),
         upp_ci = qnorm(0.975, mean = pred, sd = sqrt(var_mu))) %>%
  left_join(df_real, by = c("group", "X1")) %>%
  mutate(group = paste0("ID: ", group),
         X1 = round(X1, 1),
         X_int = cut(X1, -0.0001:9)) %>%
  group_by(group, X_int) %>%
  arrange(group, X1) %>%
  mutate(id = 1:n()) %>%
  mutate(
    low_ci = ifelse(id != 14, NA, low_ci),
    upp_ci = ifelse(id != 14, NA, upp_ci),
  )


rss_calc <- test %>% 
  filter(group %in% c(308, 309, 351)) 

pred_test_filter <- predict_hbm(model = plot_model, newdata = rss_calc, 
                         formula = formula, 
                         group_variable = group_variable)
pred_lme <- predict(lm3_m0_normal, rss_calc)
rss_hbart <- round(rss(pred_test_filter$pred, rss_calc$y), 2)
rss_lme <- round(rss(pred_lme, rss_calc$y), 2)

pred_test_final %>% 
  ggplot(aes(x = X1, y = pred)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill = "#75E6DA", alpha = 0.3) + 
  geom_point(aes(colour = "#F96209"), size = 0.5) +
  geom_point(colour = "#F96209", size = 0.5) +
  geom_point(aes(x = X1, y = y, colour = "#75E6DA"), size = 0.75) + 
  geom_point(aes(x = X1, y = y), colour = "#75E6DA", size = 2) + 
  geom_line(aes(x = X1, y = fit), colour = "#75E6DA") + 
  #geom_ribbon(aes(ymin = lwr_gam, ymax = fit_gam), fill = "lightblue") + 
  geom_errorbar(aes(ymin = low_ci, ymax = upp_ci),
                 position = position_dodge(width = 0.2),
                width = 0.5, colour = '#F96209') +
  #geom_smooth(colour = "grey80",  alpha = 0.2) +
  facet_wrap(~group, ncol = 3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(y = "Response time (ms)", 
       x = 'Covariate: days of sleep deprivation', 
  title = paste0("Residual sum of squares\nHBART: ", rss_hbart, ", LME: ", rss_lme)
  ) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_linedraw(20) +
  scale_colour_manual(
    name="Source:",
    values=c(Data="#75E6DA", 
             `HBART Prediction`="#F96209"), 
    guide = guide_legend(override.aes = list(
      size = c(3, 3), shape = c(16, 16)))) + 
  theme(panel.spacing.x = unit(0.5, "lines"), 
        legend.position = "bottom")

# ggsave(file = "results/reports/february/predictions_plot.png",
#        width = 12, height = 6)

ggsave(file = "paper/predictions_plot.png",
       width = 8, height = 5)

#-------------------------------------------------------------------------------
pred_test <- predict_hbm(model = plot_model_missing, newdata = new_test, 
                         formula = formula, 
                         group_variable = group_variable) %>% 
  filter(group %in% c(308, 309, 351)) 
new_test <- new_test %>% 
  filter(group %in% c(308, 309, 351)) 
# LME without a few groups is just lm? 
lm_missing <- lm(y ~ X1, data = train_2)

pred_test$fit <- predict(lm3_m0_missing, new_test, re.form=NA)
pred_lm <- predict(lm_missing, new_test, interval = "predict")
pred_test$upr <- pred_lm[, 3]
pred_test$lwr <- pred_lm[, 2]
#pred_test$upp <- pred_test$fit + 1.96 * 

#summary(lm3_m0_missing)



pred_test_final <- pred_test %>%
  mutate(low_ci = qnorm(0.025, mean = pred, sd = sqrt(var_mu)),
         upp_ci = qnorm(0.975, mean = pred, sd = sqrt(var_mu))) %>%
  left_join(df_real, by = c("group", "X1")) %>%
  mutate(group = paste0("ID: ", group),
         X1 = round(X1, 1),
         X_int = cut(X1, -0.0001:9)) %>%
  group_by(group, X_int) %>%
  arrange(group, X1) %>%
  mutate(id = 1:n()) %>%
  mutate(
    low_ci = ifelse(id != 14, NA, low_ci),
    upp_ci = ifelse(id != 14, NA, upp_ci),
  )


rss_calc <- test_2 %>% 
  filter(group %in% c(308, 309, 351)) 

pred_test_filter <- predict_hbm(model = plot_model, newdata = rss_calc, 
                                formula = formula, 
                                group_variable = group_variable)
pred_lme <- predict(lm3_m0_normal, rss_calc)
rss_hbart <- round(rss(pred_test_filter$pred, rss_calc$y), 2)
rss_lme <- round(rss(pred_lme, rss_calc$y), 2)


pred_test_final %>% 
  ggplot(aes(x = X1, y = pred)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill = "#75E6DA", alpha = 0.3) + 
  geom_point(aes(colour = "#F96209"), size = 0.5) +
  geom_point(colour = "#F96209", size = 0.5) +
  geom_point(aes(x = X1, y = y, colour = "#75E6DA"), size = 0.75) + 
  geom_point(aes(x = X1, y = y), colour = "#75E6DA", size = 2) + 
  geom_line(aes(x = X1, y = fit), colour = "#75E6DA") + 
  #geom_ribbon(aes(ymin = lwr_gam, ymax = fit_gam), fill = "lightblue") + 
  geom_errorbar(aes(ymin = low_ci, ymax = upp_ci),
                position = position_dodge(width = 0.2),
                width = 0.5, colour = '#F96209') +
  #geom_smooth(colour = "grey80",  alpha = 0.2) +
  facet_wrap(~group, ncol = 3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(y = "Response time (ms)", 
       x = 'Covariate: days of sleep deprivation', 
       title = paste0("Residual sum of squares\nHBART: ", rss_hbart, ", LME: ", rss_lme)
       ) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_linedraw(20) +
  scale_colour_manual(
    name="Source:",
    values=c(Data="#75E6DA", 
             `HBART Prediction`="#F96209"), 
    guide = guide_legend(override.aes = list(
      size = c(3, 3), shape = c(16, 16)))) + 
  theme(panel.spacing.x = unit(0.5, "lines"), 
        legend.position = "bottom")

ggsave(file = "paper/predictions_plot_missing.png",
       width = 8, height = 5)

# pred_test <- pred_test %>% 
#   mutate(low_ci = qnorm(0.025, mean = pred, sd = sqrt(var_mu)), 
#          upp_ci = qnorm(0.975, mean = pred, sd = sqrt(var_mu))) %>% 
#   left_join(df_real, by = c("group", "X1")) %>% 
#   mutate(group = paste0("ID: ", group), 
#          X1 = round(X1, 1), 
#          X_int = cut(X1, -0.001:9)) %>% 
#   group_by(group, X_int) %>% 
#   mutate(id = 1:n()) %>% 
#   mutate(
#     low_ci = ifelse(id != 3, NA, low_ci),
#     upp_ci = ifelse(id != 3, NA, upp_ci),
#   )
# 
# pred_test %>% 
#   ggplot(aes(x = X1, y = pred)) +
#   geom_point(aes(colour = "#F96209"), size = 0.5) +
#   geom_point(colour = "#F96209", size = 0.5) +
#   geom_point(aes(x = X1, y = y, colour = "#75E6DA"), size = 0.75) + 
#   geom_point(aes(x = X1, y = y), colour = "#75E6DA", size = 1.5) + 
#   geom_errorbar(aes(ymin = low_ci, ymax = upp_ci),
#                 position = position_dodge(width = 0.2),
#                 width = 0.5, colour = '#F96209') +
#   geom_smooth(colour = "grey80",  alpha = 0.2) +
#   facet_wrap(~group, ncol = 6) +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
#   labs(y = "Response time", 
#        x = 'Covariate: days of sleep deprivation', 
#        title = "Predictions for patients in the sleep study") + 
#   ylim(min(pred_test$pred) - 0.5, max(pred_test$pred) + 0.7) +
#   theme_linedraw(16) +
#   scale_colour_manual(
#     name="Source:",
#     values=c(Data="#75E6DA", 
#              `HBM Prediction`="#F96209"), 
#     guide = guide_legend(override.aes = list(
#       size = c(3, 3), shape = c(16, 16)))) + 
#   theme(panel.spacing.x = unit(0.5, "lines"), 
#         legend.position = "bottom")
# 
# ggsave(file = "results/reports/february/predictions_plot_missing_id.png",
#        width = 13, height = 8.5)

#-------------------------------------------------------------------------------

# split_g <- pred_test %>% split(.$group)
# split_g <- list(split_g$`308`, split_g$`309`, split_g$`351`)
# 
# pred_split <- new_test %>% 
#   filter(group %in% c(308, 309, 351)) %>% 
#   split(.$group)
# pred_split <- list(pred_split$`308`, pred_split$`309`, pred_split$`351`)
# 
# fit_gam <- function(train_gam, pred_gam, group){
#   m <- gam(pred ~ s(X1), data = train_gam) 
#   p <- predict(m, newdata = pred_gam, se = TRUE) 
# 
#   g <- data.frame(X1 = pred_gam$X1, 
#                   fit_gam = p$fit,
#                   lwr_gam = p$fit - 5.96*p$se.fit, 
#                   upr_gam = p$fit + 5.96*p$se.fit, 
#                   group = group)
#   g
# }
# 
# groups <- c(308, 309, 351)
# gg <- pmap_dfr(list(split_g, pred_split, groups), fit_gam)
# gg$group <- as.factor(gg$group)

# pred_test2 <- pred_test %>% 
#   left_join(gg, by = c("X1", "group"))

#-------------------------------------------------------------------------------