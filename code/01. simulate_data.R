# May 2021
library(tidyverse)
alpha = 0.5; beta = 1; mu_mu = 0; k1 = 8 ; k2 = 10
# 1. Simulating one very simple tree (1 node)
set.seed(2022)
n <- 2000
m <- 5
alloc <- sample(1:m, size = n, replace = TRUE)
X1 <- runif(n)
tau <- rgamma(1, 1/alpha, beta) # 1.29
mu <- rnorm(2, mu_mu, sqrt(k2/tau))
mu

muj_1 <- rnorm(m, mu[1], sqrt(k1/tau))
muj_2 <- rnorm(m, mu[2], sqrt(k1/tau))

for(i in 1:n) {
  curr_mean <- if(X1[i] < 0.5) { muj_1
  } else { muj_2 }
  y[i] <- rnorm(1, curr_mean[alloc[i]], sd = sqrt(1/tau))
}


#-----------------------------------------------------------------
data_one_node <- data.frame(X1, y, group = alloc) %>% 
  mutate(set = ifelse(runif(n()) > 0.8, "test", "train")) %>% 
  split(.$set)
#-----------------------------------------------------------------
# 2. Simulating a more complicated tree
sim_data <- function(m = 5, n = 500, n_nodes = 3, 
                     k1 = 1, k2 = 3){
  #set.seed(2025)
  alloc <- sample(1:m, size = n, replace = TRUE)
  alpha = 0.5; beta = 1; mu_mu = 0; 
  X1 <- runif(n)
  X2 <- runif(n)
  #tau <- rgamma(1, 1/alpha, beta)
  tau <- 2
  mu <- rnorm(n_nodes, mu_mu, sqrt(k2/tau))
  muj_1 <- round(rnorm(m, mu[1], sqrt(k1/tau)), 3)
  muj_2 <- round(rnorm(m, mu[2], sqrt(k1/tau)), 3)
  muj_3 <- round(rnorm(m, mu[3], sqrt(k1/tau)), 3)
  mu_js <- list(muj_1, muj_2, muj_3)
  var_values_save <- vector()
  sim_norm <- function(group, mean_list, taub = tau){
    rnorm(n = 1, mean = mean_list[group],  
          sd = sqrt(1/taub))
  }
  
  vars <- c("X1", "X2")
  y <-  rep(NA, n)
  data <- tibble(y, X1, X2, group =  alloc)
  data$parent <- 1
  data$mean_mu <- list(muj_1)
  
  for(i in 1:(n_nodes - 1)){
    sample_vars <- sample(size = 1, x = vars)
    vars <- vars[-which(vars == sample_vars)]
    sample_parent <- sample(size = 1, 
                            x = unique(data$parent))
    
    
    var_values <- data %>% pull(sample_vars)
    bp <- median(var_values)
    np <- sample_parent+1
    
    new_break <- data %>% 
      filter(parent == sample_parent) %>% 
      mutate(
        mean_mu = ifelse(.data[[sample_vars]] <= bp, 
                         mu_js[i+1], mean_mu), 
        parent = ifelse(.data[[sample_vars]] <= bp, 
                        np, parent))
    
    data <- bind_rows(new_break, 
                      data %>% 
                        filter(parent != sample_parent))
    var_values_save <- c(bp, var_values_save)
    
  }
  
  data <- data %>% 
    mutate(y = map2_dbl(group, mean_mu, sim_norm)) %>% 
    mutate(set = ifelse(runif(n()) > 0.8, "test", "train")) %>% 
    split(.$set)
  
  return(list(y = y, 
              data = data, 
              tau = tau, 
              mu = mu, 
              var_values_save = var_values_save, 
              mujs = mu_js))
  }

# Small k1 and k2 
sim_ks_small <- sim_data(n = 2000, k1 = 1, k2 = 3)

# Small k1, bigger k2
sim_k2_big <- sim_data(n = 2000, k1 = 2.5, k2 = 10)

# Small k2, bigger k1
sim_k1_big <- sim_data(n = 2000, k1 = 15, k2 = 2)

# Bigger ks
sim_ks_big <- sim_data(n = 2000, k1 = 8, k2 = 10)

saveRDS(data_one_node, "data/data_one_node.rds")

saveRDS(list(
             sim_ks_small = sim_ks_small, 
             sim_k2_big = sim_k2_big, 
             sim_k1_big = sim_k1_big, 
             sim_ks_big = sim_ks_big), 
        "data/all_data.rds")
#-----------------------------------------------------------------
library(tidyverse)
alpha = 0.5; beta = 1; mu_mu = 0; k1 = 8 ; k2 = 10
# 1. Simulating one very simple tree (1 node)
set.seed(2022)
n <- 2000
m <- 25
alloc <- sample(1:m, size = n, replace = TRUE)
X1 <- runif(n)
tau <- rgamma(1, 1/alpha, beta) # 1.29
mu <- rnorm(2, mu_mu, sqrt(k2/tau))
mu # -1.628563  1.483883

muj_1 <- rnorm(m, mu[1], sqrt(k1/tau))
muj_2 <- rnorm(m, mu[2], sqrt(k1/tau))
sd(c(muj_1, muj_2))
sqrt(k1/tau)

qplot(alloc, y) +
  geom_boxplot(aes(group = alloc))

for(i in 1:n) {
  curr_mean <- if(X1[i] < 0.5) { muj_1
  } else { muj_2 }
  y[i] <- rnorm(1, curr_mean[alloc[i]], sd = sqrt(1/tau))
}

#-----------------------------------------------------------------
data_k1 <- data.frame(X1, y, group = alloc) 
  #mutate(set = ifelse(runif(n()) > 0.9, "test", "train")) %>% 
  #split(.$set)
#tau = 1.535821
saveRDS(data_k1, "data/data_k1.rds")
