posteriors <- function(params, tau_p = NULL, 
                       mu_j_p = NULL,
                       mu_p = NULL){
  N <- params$N
  k <- params$k
  beta <- params$beta
  alpha <- params$alpha
  group <- params$group
  mu_j <- params$mu_j
  tau <- params$tau
  tau_mu <- params$tau_mu
  mu_main <- params$mu_main
  y <- params$y
  m <-  length(mu_j)
  
  if(!is.null(tau_p)){
    #alpha_tau <- (N + m)/2 + alpha
    alpha_tau <- N/2 + alpha
    term_mu <- c()
    term_mu_j <- c()
    for(j in unique(group)){
      y_j <- y[group == j]
      term_mu[j] <- sum((y_j - mu_j[j])^2)
      term_mu_j[j] <- (mu_j[j] - mu_main)^2
    }
    
    
    beta_tau <- sum(term_mu)/2 + beta + sum(term_mu_j)/(2 *k)
    alpha_tau/beta_tau
    
    post <- dgamma(tau_p, alpha_tau, beta_tau)
  }
  
  
  if(!is.null(mu_j_p)){
    mu_j_p <-  rep(mu_j_p, length(nj))
    mean_mu <- c()
    var_mu <- c()
    
    for(j in unique(group)){
      y_bar_j <- mean(y[group == j]) 
      #mean_mu[j] <- tau * ((mu_main/k) +  y_bar_j * nj[j])/(nj[j] + 1/k)
      mean_mu[j] <- ((mu_main/k) +  y_bar_j * nj[j])/(nj[j] + 1/k)
      var_mu[j] <- (1/(nj[j] + 1/k))* tau
    }
    post <- dnorm(mu_j_p, mean = mean_mu, sd = sqrt(var_mu))
  } 
  
  if(!is.null(mu_p)){
    mean_mu <- (tau/k) * mean(mu_j) * m / (tau_mu + (tau/k)*m)
    #mean_mu <- ((1/k) * mean(mu_j))/ (tau_mu + (tau/k)*m)
    var_mu <- (tau_mu + (tau/k)*m)^(-1)
    post <- dnorm(mu_p, mean = mean_mu, sd = sqrt(var_mu))
  } 
  
  return(list(post = post))
  
}
