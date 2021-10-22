library(tidyverse)
files <- list.files("code/mixedbart/R") %>% 
  paste0("code/mixedbart/R/", .)
map(c(files[-c(7, 10, 11)]), source)
rss <- function(y, y_test){ mean((y - y_test)^2)}

data <- list.files("data")[4] %>% 
  paste0("data/", .) %>% 
  read_csv()
glimpse(data)

# View(data)
# names(data)
#------------------------------------------------------------
iter <-  2500
lim <- 1500
group_variable <- "Rep"
formula <- (compr_strength ~ mean_pd1 + mean_pd2 + mean_bd +
              quantile_1_pd1 + quantile_1_pd2 + quantile_1_bd +
              quantile_99_pd1 + quantile_99_pd2 + quantile_99_bd) %>% 
  as.formula()

alpha = 0.5; beta = 1; mu_mu = 0;

data$compr_strength <- c(scale(data$compr_strength))

pars <- list(
  k1 = 5, k2 = 3, alpha = alpha, beta = beta, 
  mu_mu = mean(data$compr_strength)
)


m0  <- bcart(formula, 
             dataset = data, 
             iter = iter, 
             group_variable, pars, scale = FALSE)

1/m0$tau_post[5000]

View(m0$final_tree)
rss(m0$final_tree$y, m0$final_tree$mu_js_sampled) %>% 
  sqrt()


library(rpart)
m1 <- rpart(formula, data)
rpart.plot::rpart.plot(m1)
rss(data$compr_strength, predict(m1)) %>% 
  sqrt()

library(ranger)
m2 <- ranger(formula, data)
rss(data$compr_strength, m2$predictions) %>% 
  sqrt()
