#' @name prune_tree
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Prunes the current tree
#' @param current_tree The current tree
#' @param drawn_node The node to prune
#' @param variable_in_question Variable that was split in the node that will
# be now pruned
#' @param nodes_to_prune The two nodes (left and right) to be pruned
#' @return The new tree
#' @example

prune_tree <- function(current_tree, drawn_node, variable_in_question){
  nodes_to_prune <- stringr::str_remove(drawn_node, '( right| left)$')
  # Updating the node that will suffer the prune (returning
  # to the point where it was before the growing)
  new_node <- stringr::str_remove(nodes_to_prune, '( X[0-9])$')
  
  current_tree %>%
    dplyr::mutate(
      # Calculating the reduction in the depth of the pruned
      # node
      d_reduction =
        ifelse(
          stringr::str_detect(node, nodes_to_prune),
          stringr::str_match(node,
                             paste0('(?<=', nodes_to_prune, '\\s).*')) %>%
            stringr::str_count(pattern = "X[0-9]") + 1, 0),
      
      # Reducing the depth of the nodes
      d = ifelse(
        stringr::str_detect(node, nodes_to_prune), d - d_reduction, d),
      
      # Changing the node for the new_node (with the pruning)
      temp_node = ifelse(
        stringr::str_detect(node, nodes_to_prune), paste(new_node), node),
      
      # Removing from the new_node the last occurrence of the name of
      # a variable + the side of the node, to update the parent
      parent = ifelse(
        stringr::str_detect(node, nodes_to_prune),
        stringr::str_remove(temp_node,
                            paste0('( X[0-9] right| X[0-9] left)$')), parent),
      
      # Updating the node and node index
      node = temp_node,
      node_index = as.numeric(as.factor(node))) %>%
    # Discarding temporary variables
    dplyr::select(-temp_node, -d_reduction)
}

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
                                   tree, current_node, p_split,
                                   var_in_prune, i, 
                                   p_grow){
  #p_grow = 0.5
  p_prune = 1 - p_grow
  # Number of available final nodes to prune -------
  b <-  dplyr::n_distinct(old_tree$node_index)
  
  # Number of internal nodes  -----------------------
  w_2 <- i
  # Probability of pruning -------------------------
  p_t_to_tstar <- p_prune/w_2

  # Probability of splitting a variable ------------
  p_adj <- p_split
  
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

  # filtered_tree <- old_tree %>%
  #   dplyr::filter(stringr::str_detect(node, nodes_to_prune)) %>%
  #   dplyr::mutate(node =
  #                   ifelse(
  #                     stringr::str_detect(node, paste0(nodes_to_prune, " left")),
  #                     paste0(nodes_to_prune, " left"),
  #                     paste0(nodes_to_prune, " right")
  #                   ))
  # The first node is on the left, the second is on the right,
  # meaning that the left node has the smaller index --------

  
  #cond_parent <- cond_calculation(data_cond = tree, pars = pars)
  #tree_b <- tree %>% filter()
  nam <- unique(tree$node)
  cond_new <- 0
  
  for(name in  nam){
    data_set <- tree %>% 
      dplyr::filter(node == name)
    marg <- cond_calculation(data_cond = data_set, pars = pars)
    cond_new <- marg + cond_new
  }
  
  nam <- unique(old_tree$node)
  cond_parent <- 0
  
  for(name in  nam){
    data_set <- old_tree %>% 
      dplyr::filter(node == name)
    marg <- cond_calculation(data_cond = data_set, pars = pars)
    cond_parent <- marg + cond_parent
  }
  
  
  lk_ratio <-  cond_new - cond_parent
  
  
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
                                  var_in_prune, p_split){
  
  # Finding the probability of selecting one
  # available predictor -------------------------------------
  p_adj <- 1/p_split
  
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
                        p_split, var_in_prune, nodes_to_prune,
                        i, mu_res_cond, p_grow){
  # All ratios:
  trans <- transition_ratio_prune(old_tree, tree, current_node,
                                  var_in_prune = var_in_prune,
                                  p_split = p_split, i = i, 
                                  p_grow = p_grow)
  
  lk <- lk_ratio_prune(old_tree, tree, current_node, pars = pars,
                       nodes_to_prune = nodes_to_prune)
  struct <- structure_ratio_prune(old_tree, tree, current_node, var_in_prune,
                                  p_split = p_split)

  r <- min(1, exp(trans+lk+struct))
  return(r)
}


