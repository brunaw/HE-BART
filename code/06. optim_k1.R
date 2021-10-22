#------------------------------------------------------------
#
#
#------------------------------------------------------------
library(tidyverse)
library(MASS)
#------------------------------------------------------------
n <- 1000
m <- 5
set.seed(2021)
x <- factor(sample(1:m, size = n, replace = TRUE))
n_j <- table(x)
M <- model.matrix(~ x - 1)
mu_mu <- 3
alpha <- 0.5
beta <- 1
k_1 <- 5
k_2 <- 8
ones_n <- rep(1, n)
MtM <- tcrossprod(x = M) # Note tcrossprod(X) computes X%*%t(X)
tMM <- crossprod(x = M) # crossprod(X) computes t(X)%*%X
W_1 <- k_1 * MtM + k_2 * tcrossprod(x = ones_n) + diag(n) 

tau <- rgamma(1, alpha, beta)
y <- mvrnorm(1, rep(mu_mu, n), W_1/tau)


inv2 <- function(k_1, k_2, M) {
  n <- nrow(M)
  n_j <- colSums(M)
  Psi_tilde_inv <- diag(n) - M%*%diag(k_1/(1 + k_1*n_j))%*%t(M)
  k_3 <- 1/k_2 + sum(Psi_tilde_inv)
  Psi_row_sums <- rowSums(Psi_tilde_inv)
  Psi_col_sums <- colSums(Psi_tilde_inv)
  return(Psi_tilde_inv - tcrossprod(Psi_row_sums, Psi_col_sums) / k_3)
}

marg2 <- function(k_1, 
                  k_2f = k_2, yf = y, 
                  Mf = M, mu_muf = mu_mu, alphaf = alpha,
                  betaf = beta){
  n <- nrow(Mf)
  W_0 <- rep(mu_muf, n)
  ymW_0 <- yf - W_0
  p_weibull <- dweibull(k_1, shape = 3, 3, log = TRUE)
  p_total <- -(
    p_weibull 
    -(n/2)*log(2*pi) - 0.5 * det2(k_1, k_2f, Mf) 
    + lgamma(n/2 + alphaf) - (n/2 + alphaf)*
      log(0.5 * t(ymW_0)%*%inv2(k_1, k_2f, Mf)%*%ymW_0 + betaf))
  return(p_total)
}
#marg2(y = y, k_1 = k_1, k_2 = k_2, M = M, mu_mu = mu_mu, alpha = alpha, beta = beta)

marg2(k_1 = 3)
marg2(k_1 = 5)
tryy <- optim(3, marg2, lower = 0, upper = 50,
              #method = "L-BFGS-B", 
              control = list(factr = 1e-10, maxit = 50))
tryy
marg2(k_1 = tryy$par[1])
marg2(k_1 = 3)
exp(3693 - 3693.504)

#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
# Optim gets resuls similar to Andrew's ones

#-----------------------------------------------------------------
group_variable = "group"
formula <- y ~ X1
iter <-  500
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

tree <- m0_2node$final_tree
k1 <- pars$k1
k2 <- pars$k2

# ------------------- this might be for each node?
# error_y <- tree %>% 
#   dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
#   dplyr::summarise(sum_errors_y = sum(err_y))


cond_calculation <- function(k1, data_cond = tree, 
                             k2 = 10, 
                             tau_post = 1.35, M = 5){
  group = data_cond$group
  n_nodes <- n_distinct(data_cond$node)
  mu_post <- unique(data_cond$mu_sampled)
  if(length(mu_post) == 1){
    mu_post <- rep(mu_post, n_nodes)
  }
  njs <-  data_cond %>%
    dplyr::group_by(node, group) %>%
    dplyr::count()
  split_nodes <- njs %>% split(.$node)
  
  n_j_means <- data_cond %>%
    group_by(node, group) %>%
    summarise(mean_y_j = mean(y)) %>%
    arrange(group) %>%
    split(.$node)
  
  mu_js_post <- list()
  nodes_unique <- unique(njs$node)
  
  for(m in 1:n_nodes){
    mean_mu <- c()
    var_mu <- c()
    
    nj_node <- split_nodes[[nodes_unique[m]]] %>% dplyr::pull(n)
    y_bar_node <-  n_j_means[[nodes_unique[m]]] %>% dplyr::pull(mean_y_j)
    
    for(j in sort(unique(group))){
      y_bar_j <- y_bar_node[j]
      mean_mu[j] <- ((mu_post[m]/k1) +  y_bar_j * nj_node[j])/(nj_node[j] + 1/k1)
      var_mu[j] <- (tau_post*(nj_node[j] + 1/k1))^(-1)
    }
    mu_js_post[[m]] <- stats::rnorm(M, mean = mean_mu, sd = sqrt(var_mu))
  }
  
  error <- expand_grid(group = 1:5,
                           node = unique(data_cond$node)) %>%
    arrange(node, group) %>%
    mutate(muj = c(mu_js_post[[1]], mu_js_post[[2]])) %>%
    left_join(data_cond %>% select(y, node, group),
              by = c("node", "group")) %>%
    dplyr::mutate(err_y = (y - muj)^2) %>%
    dplyr::summarise(sum_errors_y = sum(err_y)) %>%
    pull(sum_errors_y)
  
  
  #y <- data_cond$y
  M_mat <- stats::model.matrix(y ~ factor(group) - 1, data_cond)
  N <- nrow(data_cond)
  psi <- (k1 * M_mat %*% t(M_mat)) + diag(N)
  #W_0 <- rep(mu_mu, N)
  W_1 <- (k2 * diag(x = 1, nrow = N, ncol = 1) %*%
            t(diag(x = 1, nrow = N, ncol = 1))) + psi
  inner <- error/2 + beta
  p_cond_y <- (-1/2) *
    log(det(W_1)) + (log(inner)*(-(N/2 + alpha))) +
    lgamma(N/2 + alpha)
  -p_cond_y
}



set.seed(2021)
seq_k1 <- abs(seq(0.0005, 15, by = 0.5) + rnorm(15))
llks <- map_dbl(seq_k1, cond_calculation)
# For this data the small values don't get picked?

data.frame(llks, llks) %>% 
  ggplot(aes(seq_k1, llks)) +
  geom_point()

cond_calculation(seq_k1[1])
cond_calculation(25)

