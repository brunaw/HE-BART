#' @name transition_ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Transition ratio for grow.
#' @description Transition probability of going to the candidate tree,
#' given that the action step was a grow.
#' @param tree The current tree.
#' @param current_node The current node.
#' @param p The number of predictors still available.
#' @param current_selec_var The variable selected for the split.
#' @param results_f The current results to find the number of second
#' generation internal nodes.
#' @return The transition ratio
#' @details When transitioning to a grow, we need the probabilities of:
#'  1. Growing the tree
#'  2. Growing from the specific node, that has to consider the
#'  available predictors and available values for that grow.
#' When transitioning from a grow back to the split,
#' we need the probabilities of:
#' 1. Pruning the tree
#' 2. Selecting  node to prune
#' @example
#' transition_ratio_grow(tree, current_node)

transition_ratio_grow <- function(tree, current_node, p,
                                  current_selec_var, results_f, 
                                  p_grow){

  #p_grow = 0.5
  p_prune = 1 - p_grow
  # Number of available final nodes to break -------
  b <-  tree %>% dplyr::distinct(node_index) %>% nrow()

  # Probability of splitting (n variables) ------------
  p_adj <-  1/p

  # Available values to split given the variable selected for the
  # split and the node selected ----------------------
  n_j_adj <-  tree %>%
    dplyr::filter(parent == current_node) %>%
    dplyr::distinct(!!rlang::sym(current_selec_var)) %>% nrow()

  
  # Probability of the transition -------------------
  p_t_to_tstar <- p_grow * (1/b) * (p_adj) * (1/n_j_adj)


  #  Probability of transitioning from the new tree
  # back to the original -------------------------
  w_2 <-  nrow(results_f)
  p_tstar_to_t <- p_prune/w_2

  trans_ratio <- log(p_t_to_tstar/p_tstar_to_t)
  return(trans_ratio)
}

#' @name lk_ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Likelihood ratio for grow.
#' @description Likelihood ratio of the candidate tree and the previous
#' tree, given that the action step was a grow.
#' @param tree The current tree.
#' @param current_node The current pruned node.
#' @param pars The full list of parameters.
#' @return The likelihood ratio.
#' @details The likelihood ratio for growth needs the
#' joint distribution of each region (of each node), given the
#' parameter sigma^2 (here, mu is considered to be integrated out)
#' @example
#' lk_ratio_grow(tree, current_node, sigma_2_y, sigma_2_mu)
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


# Full parent conditional -----------------------------
cond_calculation <- function(data_cond, error = 0, pars){
  # y <- data_cond$y
  # M_mat <- stats::model.matrix(y ~ factor(group) - 1, data_cond)
  # N <- nrow(data_cond)
  # psi <- (k1 * M_mat %*% t(M_mat)) + diag(N)
  # W_0 <- rep(mu_mu, N)
  # W_1 <- (k2 * diag(x = 1, nrow = N, ncol = 1) %*%
  #           t(diag(x = 1, nrow = N, ncol = 1))) + psi
  # 
  # inner <- error/2 + beta
  # p_cond_y <- (-1/2) *
  #   log(det(W_1)) + (log(inner)*(-(N/2 + alpha))) +
  #   lgamma(N/2 + alpha)
  
  # error <- data_cond %>% 
  #   mutate(error = (y - mu_js_sampled)^2) %>% 
  #   summarise(error = sum(error)) %>% 
  #   pull(error)
  
  M <- model.matrix(~ factor(data_cond$group) - 1)
  y <- data_cond$y
  k_1 = pars$k1 
  k_2 = pars$k2 
  mu_mu = pars$mu_mu 
  alpha = pars$alpha
  beta = pars$beta
  
  n <- nrow(M)
  #n <- 1
  W_0 <- rep(mu_mu, n)
  ymW_0 <- y - W_0
  
  term_1 <- -(n/2)*log(2*pi)
  term_2 <- - 0.5 * det2(k_1_d = k_1, k_2_d = k_2, M_d = M)
  term_3 <- lgamma(n/2 + alpha)
  
  term_4 <- - (n/2 + alpha)*log(0.5 * t(ymW_0)%*%inv2(k_1, k_2, M)%*%ymW_0 + beta)
  #term_4 <- - (n/2 + alpha)*log(0.5 * error + beta)
  #p_cond_y <- term_1 + term_3 + term_2 + term_4
  p_cond_y <- term_1 + term_4 + term_2 + term_3
  
}



lk_ratio_grow <- function(tree, current_node, pars){
  
  filtered_tree <- tree %>% 
    dplyr::filter(parent == current_node)

  # The first node is on the left, the second is on the right,
  # meaning that the left node has the smaller index --------

  # # Counting how many observations are in each
  # # region (left and right) ---------------------------------
  # nl_nr <-  filtered_tree %>%
  #   dplyr::count(node_index) %>%
  #   dplyr::arrange(n) %>%
  #   dplyr::pull(n)
  #
  # # Counting how many observations are in each group and node
  # # region (left and right) ---------------------------------
  # n_group_node <-  filtered_tree %>%
  #   dplyr::count(group, node_index)
  #
  # # Calculating the sums of y in each region ----------------
  # sums_nodes <- filtered_tree %>%
  #   dplyr::group_by(node_index) %>%
  #   dplyr::summarise(sum = sum(y)) %>%
  #   dplyr::arrange(node_index) %>%
  #   dplyr::pull(sum)


  # beta <- pars$beta
  # alpha <- pars$alpha
  # mu_mu <- pars$mu_mu
  # k1 <- pars$k1
  # k2 <- pars$k2
  
  #---------------------------------------------
  # Conditional split by node
  # split_node <- filtered_tree %>%
  #   split(.$criteria)
  # 
  
  # error_y <- filtered_tree %>%
  #   dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
  #   dplyr::summarise(sum_errors_y = sum(err_y))

  #n_nodes <- n_distinct(filtered_tree$node)
  cond_parent <- cond_calculation(data_cond = filtered_tree, pars = pars)

  nam <- unique(filtered_tree$node)
  cond_node <- 0
  
  for(i in 1:length(nam)){
    data_set <- tree %>% 
      filter(node == nam[i])
    marg <- cond_calculation(data_set, pars = pars)
    cond_node <- marg + cond_node
  }
  
  lk_ratio <-  cond_node - cond_parent
  

  # exp(lk_ratio) # exploded
  # # Calculating the equation of the lk. ratio ---------------
  # first_term <- log(sqrt(
  #   ((sigma_2_y*(sigma_2_y + sigma_2_mu * sum(nl_nr)))/
  #      ((sigma_2_y + sigma_2_mu * nl_nr[1])*
  #      (sigma_2_y + sigma_2_mu * nl_nr[2]))
  #   )
  # ))
  #
  # # Exponential part -----------------------------------------
  # first_term_exp <- sigma_2_mu/(2*sigma_2_y)
  # # Left node is in the first position of the objects
  # second_term_exp <- (sums_nodes[1]^2)/(sigma_2_y + nl_nr[1] * sigma_2_mu)
  # # Right node is in the second position of the objects
  # third_term_exp <- (sums_nodes[2]^2)/(sigma_2_y + nl_nr[2] * sigma_2_mu)
  #
  # fourth_term_exp <- (sum(sums_nodes)^2)/(sigma_2_y + sum(nl_nr) * sigma_2_mu)
  #
  #
  # # The exponential part
  # # check if the - sign is there or not
  # exp_part <- first_term_exp *
  #   (second_term_exp + third_term_exp - fourth_term_exp)
  #
  # lk_ratio <- exp_part + first_term
  return(lk_ratio)
}

#' @name structure_ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Tree structure ratio for grow
#' @description Tree structure ratio of the candidate tree and
#' the previous tree, given that the action step was a grow.
#' @param tree The current tree.
#' @param current_node The current grown node.
#' @param current_selec_var The variable selected for the split.
#' @param p The number of available predictors.
#' @return The tree structure ratio.
#' @details For the tree structure ratio of the new tree, we need
#' the probabilities of:
# 1. Splitting at node n
# 2. Splitting at the node of the left
# 3. Splitting at the node of the right
# 4. Using each rule at node n
#' @example
#' structure_ratio_grow(tree, current_node)

structure_ratio_grow <- function(tree, current_node,
                                 current_selec_var, p){

  # Finding the probability of selecting one
  # available predictor -------------------------------------
  p_adj <- 1/p

  # Counting the distinct rule options from
  # this available predictor -------------------------------
  n_j_adj <-  tree %>%
    dplyr::filter(parent == current_node) %>%
    dplyr::distinct(!!rlang::sym(current_selec_var)) %>% nrow()

  # Calculating the probability of the chosen rule --------
  p_rule <- p_adj * (1/n_j_adj)

  # Calculating the probability of split
  terminal_nodes <- tree %>% dplyr::distinct(node_index) %>% nrow()
  
  # !!! this is a hack, that should be fixed to no split 
  if(terminal_nodes == 1){
    p_split = 0.5
  } else {
    p_split <- 1/terminal_nodes
  } 
  
  p_t_star <- ((1-p_split)^2)*p_split*p_rule

  p_t <- (1 - p_split)

  st_ratio <- log(p_t_star/p_t)

  return(st_ratio)
}

#' @name ratio_grow
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Final ratio for a growth step.
#' @description The final ratio is to be used as the acceptance
#' criteria in the MCMC of the b-cart model.
#' @param tree The current tree.
#' @param current_node The current grown node.
#' @param pars The full list of parameters.
#' @param p The number of available predictors
#' @param current_selec_var The variable selected for the split.
#' @param results_f The current results to find the number of second
#' generation internal nodes.
#' @return The final ratio for the candidate tree.
#' @example
#' ratio_grow(tree, current_node, sigma_2_mu, sigma_2)

ratio_grow <- function(tree, current_node,
                       pars,
                       p, current_selec_var,
                       results_f,
                       p_grow){
  # All ratios:
  trans <- transition_ratio_grow(tree, current_node,
                                 current_selec_var = current_selec_var,
                                 p = p, results_f = results_f,
                                 p_grow = p_grow)
  # This is too big 
  lk <- lk_ratio_grow(tree, current_node, pars)
  # -13.04714

  struct <- structure_ratio_grow(tree, current_node,
                                 current_selec_var = current_selec_var,  p = p)

  r <- min(1, exp(trans + lk + struct))
  return(r)
}
