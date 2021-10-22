#------------------------------------------------------
# This code creates the new functions used
# to simulate data (only returns tau and the data)
# August,  2021
#------------------------------------------------------
sim_y <- function(M, mu_mu = 0, taud = tau, k1d = k1, k2d = k2){
  n <- nrow(M)
  ones_n <- rep(1, n)
  MtM <- tcrossprod(x = M) # Note tcrossprod(X) computes X%*%t(X)
  tMM <- crossprod(x = M) # crossprod(X) computes t(X)%*%X
  W_1 <- k1d * MtM + k2d * tcrossprod(x = ones_n) + diag(n) 
  y <- mvrnorm(1, rep(mu_mu, n), W_1/taud)
  y
}
sim_a <- function(m, n = 1000){
  alpha = 0.5; beta = 1; mu_mu = 0; k1 = 8 ; k2 = 10
  tau <- rgamma(1, 1/alpha, beta)
  alloc <- factor(sample(1:m, size = n, replace = TRUE))
  X1 <- runif(n)
  group_left  <- factor(alloc[X1 < 0.5])
  group_right <- factor(alloc[!(X1 < 0.5)])
  M_1_left    <- model.matrix(~ group_left - 1)
  M_1_right   <- model.matrix(~ group_right - 1)
  
  y_left <- sim_y(M_1_left, taud = tau, k1d = k1, k2d = k2)
  y_right <- sim_y(M_1_right, taud = tau, k1d = k1, k2d = k2)
  
  data <- data.frame(y = y_left, group = group_left, X1 = X1[X1 < 0.5]) %>% 
    bind_rows(data.frame(y = y_right, group = group_right, X1 = X1[!(X1 < 0.5)]))
  return(list(data = data, tau = tau)) 
}

sim_b <- function(m, n = 1000, k1 = 8, k2 = 5){
  alpha = 0.5; beta = 1; mu_mu = 0; 
  alloc <- sample(1:m, size = n, replace = TRUE)
  X1 <- runif(n)
  tau <- rgamma(1, 1/alpha, beta) # 1.29
  mu <- rnorm(2, mu_mu, sqrt(k2/tau))
  
  muj_1 <- rnorm(m, mu[1], sqrt(k1/tau))
  muj_2 <- rnorm(m, mu[2], sqrt(k1/tau))
  y <- y2 <- c()
  for(i in 1:n) {
    curr_mean <- if(X1[i] < 0.5) { muj_1
    } else { muj_2 }
    y[i] <- rnorm(1, curr_mean[alloc[i]], sd = sqrt(1/tau))
  }
  
  for(i in 1:n) {
    curr_mean <- muj_1
    y2[i] <- rnorm(1, curr_mean[alloc[i]], sd = sqrt(1/tau))
  }
  
  data <- data.frame(X1, y, group = alloc) 
  return(list(data = data, tau = tau, y2 = y2)) 
}


sim_stump <- function(m, n = 1000, k1 = 8 ){
  alpha = 0.5; beta = 1; mu_mu = 0; k2 = 10
  alloc <- sample(1:m, size = n, replace = TRUE)
  tau <- rgamma(1, 1/alpha, beta) # 1.29
  mu <- rnorm(2, mu_mu, sqrt(k2/tau))
  X1 <- runif(n)
  muj_1 <- rnorm(m, mu[1], sqrt(k1/tau))
  
  curr_mean <- muj_1  
  for(i in 1:n) {
    y[i] <- rnorm(1, curr_mean[alloc[i]], sd = sqrt(1/tau))
  }
  
  data <- data.frame(y, X1, group = alloc) 
  return(list(data = data, tau = tau)) 
}
