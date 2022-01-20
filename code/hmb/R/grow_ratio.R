#' @name grow_tree
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Grows the current tree
#' @param current_tree The current tree
#' @param selec_var The node to split on node
#' @param drawn_node The node to split on
#' @param rule The variable value to split on
#' @return The new tree
#' @example
#' grow_tree()

grow_tree <- function(current_tree, selec_var, drawn_node, rule){
  current_tree %>%
    dplyr::mutate(
      # Increasing the depth of the node 
      d =  ifelse(node == drawn_node, d + 1, d),
      # Updating the parent of the split node
      parent = ifelse(node == drawn_node, drawn_node, parent),
      # Changing the node "side" of each observation: left and right
      criteria = ifelse(
        node == drawn_node,
        ifelse(!!rlang::sym(selec_var) > rule, "left", "right"), "no split"),
      
      # Updating the node with the new split
      node = ifelse(node == drawn_node,
                    ifelse(!!rlang::sym(selec_var) > rule,
                           paste(node, selec_var, "left"),
                           paste(node, selec_var, "right")), node),
      # Updating the node index
      node_index =  as.numeric(as.factor(node)))
}

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
                                  p_grow, i){

  #p_grow = 0.5
  p_grow = 0.95
  p_prune = 1 - p_grow
  # Number of available final nodes to split on -------
  b <-  length(unique(tree$node_index))

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
  #w_2 <-  nrow(results_f)
  w_2 <- i
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
cond_calculation <- function(data_cond, pars){
  
  M <- model.matrix(~ factor(data_cond$group) - 1)
  y <- data_cond$y
  k_1 = pars$k1 
  k_2 = pars$k2 
  mu_mu = pars$mu_mu 
  alpha = pars$alpha
  beta = pars$beta
  
  n <- nrow(M)
  W_0 <- rep(mu_mu, n)
  ymW_0 <- y - W_0
  
  term_1 <- -(n/2)*log(2*pi)
  term_2 <- - 0.5 * det2(k_1_d = k_1, k_2_d = k_2, M_d = M)
  term_3 <- lgamma(n/2 + alpha)
  term_4 <- - (n/2 + alpha)*log(0.5 * t(ymW_0)%*%inv2(k_1, k_2, M)%*%ymW_0 + beta)
  p_cond_y <- term_1 + term_4 + term_2 + term_3
  p_cond_y
}



lk_ratio_grow <- function(tree, current_node, pars){
  
  filtered_tree <- tree %>% dplyr::filter(parent == current_node)
  cond_parent <- cond_calculation(data_cond = filtered_tree, pars = pars)

  nam <- unique(filtered_tree$node)
  cond_node <- 0
  
  for(name in  nam){
    data_set <- tree %>% 
      dplyr::filter(node == name)
    marg <- cond_calculation(data_cond = data_set, pars = pars)
    cond_node <- marg + cond_node
  }
  
  lk_ratio <-  cond_node - cond_parent
  
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
                       pars, p, current_selec_var,
                       i, p_grow){
  # All ratios:
  trans <- transition_ratio_grow(tree, current_node,
                                 current_selec_var = current_selec_var,
                                 p = p, i = i,
                                 p_grow = p_grow)
  
  lk <- lk_ratio_grow(tree, current_node, pars)

  struct <- structure_ratio_grow(tree, current_node,
                                 current_selec_var = current_selec_var,  p = p)

  r <- min(1, exp(trans + lk + struct))
  return(r)
}
