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

sim_b <- function(m = 10, n = 1000, k1 = 8, k2 = 5){
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

sim_bart <- function(n, k1 = 8, m = 10, k2 = 5){
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
  return(list(data = data, tau = tau, y2 = y2, 
              mu = mu, muj_1 = muj_1, muj_2 = muj_2)) 
}



#------------------------------------------------------

sim_friedman_bart = function(n, scale_err = 2, j = 10) {
  
  alloc <- sample(1:j, size = n, replace = TRUE)
  muj1 <- rnorm(j, 0, 10)
  muj2 <- rnorm(j, 0, 2)
  
  # Simulate some data using a multivariate version of Friedman
  # y = 10sin(πx1x2)+20(x3−0.5)2+10x4+5x5+ε
  
  X = matrix(runif(n*(5+p)), nrow = n, ncol = 5)
  
  pars = c(1, 1, 1, 1)
  
  mean = pars[1] * sin(pi*X[,1]*X[,2]) + pars[2] * (X[,3]-0.5)^2 + 
    pars[3] * X[,4] + pars[4] * X[,5]
  
  y = rnorm(n, mean, scale_err) 
  
  ns <- table((X[,1] < 0.5) * 1)
  to_add1 <- c()
  to_add2 <- c()
  for(i in 1:ns[1]){
    to_add1[i] <- rnorm(1, muj1[alloc[i]], sd = 1)
  }

  for(i in 1:ns[2]){
    to_add2[i] <- rnorm(1, muj2[alloc[i]], sd = 1)
  }
  

  y[X[,1] < 0.5] <- y[X[,1] < 0.5] + to_add2
  y[!(X[,1] < 0.5)] <- y[!(X[,1] < 0.5)] + to_add1
  
  df <- data.frame(y, X, group = alloc)
  return(df)
}
