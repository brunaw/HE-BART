#' sample_parameters_bart
#' @name sample_parameters_bart
#' @rdname sample_parameters_bart
#' @param P The number of trees
#' @param i The iteration 
#' @param N The size of the data
#' @param beta A parameter from the tau prior
#' @param alpha A parameter from the tau prior
#' @param current_tree_post The current trees object
#' @param k1 The current value of the k1 parameter
#' @param k2 The current value of the k2 parameter
#' @param tau_post The current value of the posterior tau
#' @param J The number of groups
#' @param mu_post The current set of mu posterior values
#' @export
#' 
sample_parameters_bart <- function(
  type = "tau", P = P, i = i, N = NULL, beta = NULL, alpha = NULL,
  current_tree_post = NULL, k1 = NULL,  k2 = NULL, tau_post, J = NULL, 
  mu_post = NULL){
  
  
  if(type == "tau"){
    # current_tree_tau <- current_tree_post %>% 
    #   dplyr::mutate(current_tree = map(tree_data, ~{ slice(.x, (i+1))})) %>% 
    #   dplyr::select(tree_index, current_tree) %>% 
    #   tidyr::unnest(current_tree) %>% 
    #   dplyr::select(tree_index, est_tree) %>% 
    #   tidyr::unnest(est_tree) %>% 
    #   dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>% 
    #   ungroup()
    
    res <- current_tree_post %>% 
      dplyr::select(tree_index, y, sampled_mu_j) %>% 
      group_by(tree_index) %>%
      mutate(rn = row_number()) %>% 
      pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>% 
      dplyr::select(-rn) %>% 
      mutate(pred = rowSums(.[2:ncol(.)]))
    
    res_muj <- current_tree_post %>% 
      group_by(tree_index, group, node) %>% 
      summarise(sampled_mu_j = unique(sampled_mu_j),
                sampled_mu = unique(sampled_mu), 
                res = (sampled_mu_j - sampled_mu)**2, 
                .groups = 'drop')
    
    res_mu <- current_tree_post %>% 
      group_by(tree_index, node) %>% 
      summarise(sampled_mu = first(sampled_mu), 
                res = (sampled_mu)**2,
                .groups = 'drop')
    
    fits <- res$pred         
    res_sum <- sum((res$y - fits)^2)
    res_muj_sum <- sum(res_muj$res)
    res_mu_sum <- sum(res_mu$res)
    
    alpha_tau <- (N + nrow(res_muj) + nrow(res_mu))/2 + alpha 
    
    beta_tau <- (res_sum + res_muj_sum *(P/k1) + res_mu_sum*(P/k2))/2 + beta
    sampled <- stats::rgamma(1, alpha_tau, beta_tau)
    
    
  } else if(type == "mu"){
    splits <- split(current_tree_post, current_tree_post$node) 
    calc_mu <-  purrr::map_dfr(
      splits, calculate_parts_mu, k1l = k1, k2l = k2, Pp = P) %>% 
      dplyr::mutate(node = names(splits), mean = p1_mu/p2_mu)
    
    calc_mu$sampled_mu <- rnorm(
      length(splits),
      mean = calc_mu$mean,
      sd = sqrt(1 /(tau_post * calc_mu$p2_mu))
    )
    sampled <- dplyr::select(calc_mu, node, sampled_mu)
    
    
  } else if(type == "muj"){
    calc_mu_j <- current_tree_post %>% 
      dplyr::select(node, group, res) %>% 
      dplyr::group_by(node, group) %>% 
      dplyr::summarise(mean_res = sum(res), n = n(), .groups = 'drop') %>% 
      dplyr::left_join(mu_post, by = "node") %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(p1_muj = (P * sampled_mu)/k1 + mean_res, 
                    p2_muj = (n + P/k1), 
                    mean_samp = p1_muj/p2_muj)
    
    calc_mu_j$sampled_mu_j <- rnorm(nrow(calc_mu_j),
                             mean = calc_mu_j$mean_samp, 
                             sd = sqrt(1 / (tau_post * calc_mu_j$p2_muj))
    )
    sampled <- dplyr::select(calc_mu_j, node, group, sampled_mu_j)
  }
  
  
  return(sampled)
}


#' calculate_parts_mu
#' A function to calculate the terms for the mu posterior 
#' @name calculate_parts_mu
#' @rdname split The splitted data 
#' @param k1l The current value of the k1 parameter
#' @param k2l The current value of the k2 parameter
#' @param P The number of trees

calculate_parts_mu <- function(split, k1l = k1, k2l = k2, Pp = P){
  res = split$res
  n <- nrow(split)
  Mm <- model.matrix(~ factor(split$group) - 1)
  ones_n <- rep(1, n)
  MtM <- tcrossprod(x = Mm) # Note tcrossprod(X) computes X%*%t(X)
  tMM <- crossprod(x = Mm) # crossprod(X) computes t(X)%*%X
  Psi <- k1l * MtM + diag(n)
  #W_1 <- Psi + k2 * tcrossprod(x = ones_n) 
  
  p2_mu <- (t(ones_n) %*% solve(Psi, ones_n) + (Pp / k2l))
  p1_mu <- t(ones_n) %*% solve(Psi, res)
  
  return(data.frame(p2_mu = p2_mu, n = n, p1_mu = p1_mu))
  
}

#' The metropolis-hastings update for k1
#' @name MH_update_hbart
#' @rdname MH_update_hbart
#' @param P The number of trees
#' @param i The iteration 
#' @param N The size of the data
#' @param beta A parameter from the tau prior
#' @param alpha A parameter from the tau prior
#' @param current_tree_post The current trees object
#' @param k1 The current value of the k1 parameter
#' @param k2 The current value of the k2 parameter
#' @param tau_post The current value of the posterior tau
#' @param J The number of groups
#' @param mu_post The current set of mu posterior values
#' @export
#' 

MH_update_hbart <- function(current_tree_mh, i = i, 
  prior = prior_k1, min_u = min_u, max_u = max_u, pars = pars) {
  
  # current_tree_mh <- current_tree_mh %>% 
  #   dplyr::mutate(current_tree = map(tree_data, ~{ slice(.x, (i+1))})) %>% 
  #   dplyr::select(tree_index, current_tree) %>% 
  #   tidyr::unnest(current_tree) %>% 
  #   dplyr::select(tree_index, est_tree) %>% 
  #   tidyr::unnest(est_tree) %>% 
  #   dplyr::select(tree_index, y, group, sampled_mu_j, sampled_mu, node) %>% 
  #   dplyr::ungroup()
  
  new_k1 <- runif(1, min = min_u, max = max_u) # Or rnorm or whatever else you want to propose
  current <- calc_all_hbart(current_tree_mh, pars$k1, pars)
  candidate <- calc_all_hbart(current_tree_mh, new_k1, pars)
  
  if(prior){
    prior_current <- dweibull(pars$k1, shape = 7, 10, log = TRUE)
    prior_candidate <- dweibull(new_k1, shape = 7, 10, log = TRUE)
    log.alpha <- (candidate - current) + (prior_candidate - prior_current) # ll is the log likelihood, lprior the log prior
  } else {
    log.alpha <- candidate - current  
  }  
  
  accept <- log.alpha >= 0 || log.alpha >= log(runif(1))
  theta <- ifelse(accept, new_k1, pars$k1)
  return(theta)
}
