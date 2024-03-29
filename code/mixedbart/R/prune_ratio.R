#' @name transition_ratio_prune
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Transition ratio for prune.
#' @description Transition probability of going to the candidate tree,
#' given that the action step was a prune.
#' @param old_tree The previous tree.
#' @param tree The current tree.
#' @param current_node The current node.
#' @param p The number of available predictors.
#' @param var_in_prune The variable that was split in the node chosen to
#' be pruned.
#' @param results_f The current results to find the number of second
#' generation internal nodes.
#' @return The transition ratio
#' @details When transitioning to a prune, we need the probabilities of:
#' 1. Pruning the tree
#' 2. Selecting  node to prune
#' When transitioning from a prune back to the split,
#' we need the probabilities of:
#'  1. Growing the tree
#'  2. Growing from the specific pruned node, that has to consider the
#'  available predictors and available values for that grow.
#' @example
#' transition_ratio_prune(tree, current_node)

transition_ratio_prune <- function(old_tree,
                                   tree, current_node, p,
                                   var_in_prune, results_f, 
                                   p_grow){
  #p_grow = 0.5
  p_prune = 1 - p_grow
  # Number of available final nodes to prune -------
  b <-  old_tree %>% dplyr::distinct(node_index) %>% nrow()
  # Number of internal nodes  -----------------------
  w_2 <-  nrow(results_f)
  # Probability of pruning -------------------------
  p_t_to_tstar <- p_prune/w_2

  # Probability of splitting a variable ------------
  p_adj <- p

  # Available values to split ----------------------
  # Using the variable that was used in the node
  # selected for prune
  n_j_adj <-  old_tree %>%
    dplyr::filter(parent == current_node) %>%
    dplyr::distinct(!!rlang::sym(var_in_prune)) %>% nrow()

  # Probability of transitioning from the new tree
  # to the old one --------------------------------
  p_tstar_to_t <- p_grow * 1/((b-1) * p_adj * n_j_adj)

  # Transition ratio in log scale
  trans_ratio <- log(p_tstar_to_t/p_t_to_tstar)
  return(trans_ratio)
}

#' @name lk_ratio_prune
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
#' lk_ratio_prune(tree, current_node, sigma_2_y, sigma_2_mu)


lk_ratio_prune <- function(old_tree, tree, current_node, pars,
                           nodes_to_prune){

  beta <- pars$beta
  alpha <- pars$alpha
  mu_mu <- pars$mu_mu
  k1 <- pars$k1
  k2 <- pars$k2

  filtered_tree <- old_tree %>%
    dplyr::filter(stringr::str_detect(node, nodes_to_prune)) %>%
    dplyr::mutate(node =
                    ifelse(
                      stringr::str_detect(node, paste0(nodes_to_prune, " left")),
                      paste0(nodes_to_prune, " left"),
                      paste0(nodes_to_prune, " right")
                    ))

  # The first node is on the left, the second is on the right,
  # meaning that the left node has the smaller index --------


  cond_calculation <- function(data_cond, error){
    y <- data_cond$y
    M_mat <- stats::model.matrix(y ~ factor(group) - 1, data_cond)
    N <- nrow(data_cond)
    psi <- (k1 * M_mat %*% t(M_mat)) + diag(N)
    W_0 <- rep(mu_mu, N)
    W_1 <- (k2 * diag(x = 1, nrow = N, ncol = 1) %*%
              t(diag(x = 1, nrow = N, ncol = 1))) + psi
    #in_W_1 <- solve(W_1)
    #inner <- (t((y - W_0)) %*% in_W_1 %*% (y - W_0))/2 + beta
    inner <- error/2 + beta
    p_cond_y <- (-1/2) *
      log(det(W_1)) + (log(inner)*(-(N/2 + alpha))) +
      lgamma(N/2 + alpha)
    p_cond_y
  }
  
  #---------------------------------------------
  # Conditional split by node
  split_node <- filtered_tree %>%
    split(.$criteria)
  
  error_y <- filtered_tree %>% 
    dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
    dplyr::summarise(sum_errors_y = sum(err_y))
  
  cond_parent <- cond_calculation(filtered_tree, error_y$sum_errors_y)
  
  if(length(names(split_node)) == 1){
    errors_y <- split_node$`no split` %>% 
      dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
      dplyr::summarise(sum_errors_y = sum(err_y))
    
    cond_node <- cond_calculation(data_cond = split_node$`no split`,
                                  error = errors_y$sum_errors_y)
  } else {
    errors_y_left <- split_node$left %>% 
      dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
      dplyr::summarise(sum_errors_y = sum(err_y))
    
    errors_y_right <- split_node$right %>% 
      dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
      dplyr::summarise(sum_errors_y = sum(err_y))
    
    
    cond_node_left <- cond_calculation(data_cond = split_node$left,
                                       error = errors_y_left$sum_errors_y)
    cond_node_right <- cond_calculation(split_node$right,
                                        error = errors_y_right$sum_errors_y)
    cond_node <- cond_node_left + cond_node_right
  }
  
  lk_ratio <- cond_parent - cond_node

  return(lk_ratio)
}


#' @name structure_ratio_prune
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Tree structure ratio for prune.
#' @description Tree structure ratio of the candidate tree and
#' the previous tree, given that the action step was a prune.
#' @param old_tree The previous tree.
#' @param tree The current tree.
#' @param current_node The current pruned node.
#' @param var_in_prune The variable that was split in the node chosen to
#' be pruned.
#' @param p The number of available predictors.
#' @return The tree structure ratio.
#' @details For the tree structure ratio of the new pruned tree, we need
#' to calculate an inversion of the tree structure ratio for the grow.
#' We need the probabilities of:
#' 1. Splitting at node n
#' 2. Splitting at the node of the left
#' 3. Splitting at the node of the right
#' 4. Using each rule at node n
#' @example
#' structure_ratio_prune(tree, current_node)

structure_ratio_prune <- function(old_tree, tree, current_node,
                                  var_in_prune, p){

  # Finding the probability of selecting one
  # available predictor -------------------------------------
  p_adj <- 1/p

  # Counting the distinct rule options from
  # the pruned predictor ----------------------------------
  n_j_adj <-  old_tree %>%
    dplyr::filter(parent == current_node) %>%
    dplyr::distinct(!!rlang::sym(var_in_prune)) %>% nrow()

  # Calculating the probability of the chosen rule --------
  p_rule <- p_adj * (1/n_j_adj)

  # Calculating the probability of split
  terminal_nodes <- old_tree %>% dplyr::distinct(node_index) %>% nrow()
  p_split <- 1/terminal_nodes

  p_t <- ((1-p_split)^2)*p_split*p_rule

  p_t_star <- (1 - p_split)

  st_ratio <- log(p_t_star/p_t)

  return(st_ratio)
}


#' @name ratio_prune
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Final ratio for a prune step.
#' @description The final ratio is to be used as the acceptance
#' criteria in the MCMC of the b-cart model.
#' @param old_tree The previous tree.
#' @param tree The current tree.
#' @param current_node The current pruned node.
#' @param pars The full list of parameters.
#' @param p The number of available predictors
#' @param var_in_prune The variable that was split in the node chosen to
#' be pruned.
#' @param nodes_to_prune The nodes to prune.
#' @param results_f The current results to find the number of second
#' generation internal nodes.
#' @return The final ratio for the candidate tree.
#' @example
#' ratio_prune(tree, current_node, sigma_2_mu, sigma_2)

ratio_prune <- function(old_tree, tree, current_node, pars,
                        p, var_in_prune, nodes_to_prune,
                        results_f, mu_res_cond, p_grow){
  # All ratios:
  trans <- transition_ratio_prune(old_tree, tree, current_node,
                                  var_in_prune = var_in_prune,
                                  p = p, results_f = results_f, 
                                  p_grow = p_grow)
  lk <- lk_ratio_prune(old_tree, tree, current_node, pars = pars,
                       nodes_to_prune = nodes_to_prune)
  struct <- structure_ratio_prune(old_tree, tree, current_node, var_in_prune,
                                  p = p)

  r <- min(1, exp(trans+lk+struct))
  return(r)
}


