det2 <- function(k_1_d, k_2_d, M_d) {
  n <- nrow(M_d)
  n_j <- colSums(M_d)
  tMM <- crossprod(x = M_d)
  Psi_tilde_inv <- diag(n) - M_d%*%diag(k_1_d/(1 + k_1_d*n_j))%*%t(M_d)
  return(log(k_2_d) + log(1/k_2_d + sum(Psi_tilde_inv)) + sum(log(1 + k_1_d*n_j)))
}

inv2 <- function(k_1, k_2, M) {
  n <- nrow(M)
  n_j <- colSums(M)
  Psi_tilde_inv <- diag(n) - M%*%diag(k_1/(1 + k_1*n_j))%*%t(M)
  k_3 <- 1/k_2 + sum(Psi_tilde_inv)
  Psi_row_sums <- rowSums(Psi_tilde_inv)
  Psi_col_sums <- colSums(Psi_tilde_inv)
  inv_M <- Psi_tilde_inv - tcrossprod(Psi_row_sums, Psi_col_sums) / k_3
  return(inv_M)
}


marg2 <- function(y, k_1, k_2, M, mu_mu, alpha, beta) {
  n <- nrow(M)
  W_0 <- rep(mu_mu, n)
  ymW_0 <- y - W_0
  
  term_1 <- -(n/2)*log(2*pi)
  term_2 <- - 0.5 * det2(k_1_d = k_1, k_2_d = k_2, M_d = M)
  term_3 <- lgamma(n/2 + alpha)
  
  term_4 <- - (n/2 + alpha)*log(0.5 * t(ymW_0)%*%inv2(k_1, k_2, M)%*%ymW_0 + beta)
  all <- term_1 + term_3 + term_2 + term_4
  
  return(all)
}

calc_all <- function(current_tree, new_k1, pars){

  nam <- unique(current_tree$node)
  tot_lik <- 0
  for(i in 1:length(nam)){
    data_set <- current_tree %>% 
      filter(node == nam[i])
    M <- model.matrix(~ factor(data_set$group) - 1)
    marg <- marg2(y = data_set$y, 
                       k_1 = new_k1, 
                       k_2 = pars$k2, 
                       M = M, 
                       mu_mu = pars$mu_mu, 
                       alpha = pars$alpha, 
                       beta = pars$beta)
    tot_lik <- marg + tot_lik
  }
  
  return(tot_lik)
  
}


MH_update <- function(
  current_tree = current_tree, 
  prior = prior_k1, min_u = min_u, max_u = max_u, pars = pars) {
  
  new_k1 <- runif(1, min = min_u, max = max_u) # Or rnorm or whatever else you want to propose
  current <- calc_all(current_tree, pars$k1, pars)
  candidate <- calc_all(current_tree, new_k1, pars)
  
  if(prior == TRUE){
    prior_current <- dweibull(pars$k1, shape = 7, 10, log = TRUE)
    prior_candidate <- dweibull(new_k1, shape = 7, 10, log = TRUE)
    log.alpha <- candidate - current + prior_candidate - prior_current # ll is the log likelihood, lprior the log prior
  } else {
    log.alpha <- candidate - current  
  }  
  
  accept <- log.alpha >= 0 || log.alpha >= log(runif(1))
  theta <- ifelse(accept, new_k1, pars$k1)
  return(theta)
}


#' @name lk_ratio_k
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Likelihood ratio for prune.
#' @description Likelihood ratio of the candidate tree and the previous
#' tree, given that the action step was a prune.
#' @param old_tree The previous tree.
#' @param tree The current tree.
#' @param current_node The current pruned node.
#' @param pars The full list of parameters.
#' @param nodes_to_prune The nodes to prune.
#' @return The likelihood ratio.
#' @details For the likelihood ratio of the new pruned tree, we need
#' to calculate an inversion of the one in the grow version.
#' The value is based on the likelihood of all regions, given the
#' tree and the parameters, for the new and the previous trees.
#' @example
#' lk_ratio_k(tree, current_node, sigma_2_y, sigma_2_mu)


lk_ratio_k <- function(new_k1, k1, current_tree, pars, group,
                       tau_post = tau_post, M){
  beta <- pars$beta
  alpha <- pars$alpha
  mu_mu <- pars$mu_mu
  k2 <- pars$k2
  # Sample new mu ------
  n_nodes <- n_distinct(current_tree$node)
  mu_post <- unique(current_tree$mu_sampled)
  if(length(mu_post) == 1){
    mu_post <- rep(mu_post, n_nodes)
  }
  njs <-  current_tree %>%
    dplyr::group_by(node, group) %>%
    dplyr::count()
  split_nodes <- njs %>% split(.$node)

  n_j_means <- current_tree %>%
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
      mean_mu[j] <- ((mu_post[m]/new_k1) +  y_bar_j * nj_node[j])/(nj_node[j] + 1/new_k1)
      var_mu[j] <- (tau_post*(nj_node[j] + 1/new_k1))^(-1)
    }
    mu_js_post[[m]] <- stats::rnorm(M, mean = mean_mu, sd = sqrt(var_mu))
  }
  
  new_error <- expand_grid(group = 1:25,
                           node = unique(current_tree$node)) %>%
    arrange(node, group) %>%
    mutate(muj = c(mu_js_post[[1]], mu_js_post[[2]])) %>%
    left_join(current_tree %>% select(y, node, group),
              by = c("node", "group")) %>%
    dplyr::mutate(err_y = (y - muj)^2) %>%
    dplyr::summarise(sum_errors_y = sum(err_y)) %>%
    pull(sum_errors_y)

  # means <- lapply(mu_js_post, mean)
  #
  # for(m in 1:n_nodes){
  #   #means_j <- means$mean_mu_j[m]
  #   means_j <- means[[m]]
  #   mean_mu <- (1/k1) * means_j*M/(M/k1 + 1/k2)
  #
  #   # Sampling from the posterior distribution of mu ---------
  #   var_mu <- (tau_post[i]*(M/k1 + 1/k2))^(-1)
  #
  #   mu_post[m] <- stats::rnorm(1, mean = mean_mu, sd = sqrt(var_mu))
  # }
  # mu[[i + 1]] <- mu_post

  # --------------------------------------------------
  # The first node is on the left, the second is on the right,
  # meaning that the left node has the smaller index --------
  error_y <- current_tree %>%
    dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
    dplyr::summarise(sum_errors_y = sum(err_y)) %>%
    pull(sum_errors_y)

  y <- current_tree$y
  inner <- error_y/2 +  beta
  M_mat <- stats::model.matrix(y ~ factor(group) - 1, current_tree)
  N <- nrow(current_tree)
  W_0 <- rep(mu_mu, N)
  k2_mat <- (k2 * diag(x = 1, nrow = N, ncol = 1) %*%
              t(diag(x = 1, nrow = N, ncol = 1)))
  #------------------------------------------------------------------
  psi <- (k1 * M_mat %*% t(M_mat)) + diag(N)
  W_1 <- k2_mat + psi

  p_cond_current <- (-1/2) * log(det(W_1)) + (
    log(inner)*(-(N/2 + alpha)))
  
  # Candidates ------
  psi_candidate <- (new_k1 * M_mat %*% t(M_mat)) + diag(N)
  W_1_candidate <- k2_mat + psi_candidate
  inner_new <-  new_error/2 + beta

  # The determinant is much bigger
  p_cond_new <- (-1/2)*log(det(W_1_candidate)) +(
    log(inner_new)*(-(N/2 + alpha)))
  
  #p_cond_new <- log(inner_new)*(-(N/2 + alpha))

  #p_cond_current <- log(inner)*(-(N/2 + alpha))
  
  # in_W_1 <- solve(W_1_candidate)
  # inner_new <- (t((y - W_0)) %*% in_W_1 %*% (y - W_0))/2 + beta

  # p_cond_new <- (-1/2)*log(det(W_1_candidate)) +(
  #   log(inner_new)*(-(N/2 + alpha)))
  
  #sqrt(k1/tau_post)
  #p_cond_current <- (-1/2)*log(det(W_1))
  
  # +(log(inner_new)*(-(N/2 + alpha)))+
  #   lgamma(N/2 + alpha)
  # --------------------------------------------------------------------

  lk_ratio <- p_cond_new - p_cond_current
  # The Unif is just not accepting the good values? 
  r <- min(1, exp(lk_ratio))
  return(r)
}
