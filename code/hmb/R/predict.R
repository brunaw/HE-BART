# results2 <- all_results[[3]]
# # results2[3, ] <- all_results[[3]][4, ]
# # results2[4, ] <- all_results[[3]][3, ]
# results <- results2
# 
# tree <- model$tree_data[[3]] %>% tail(1)
# 
results <- m0_real$final_trees$results[2]]
tree <- m0_real$final_trees$tree_data[[2]] %>% tail(1)
table(tree$node)

calculate_preds <- function(results, tree, newdata){
  
  tree <- tree %>% 
    dplyr::select(est_tree) %>% 
    tidyr::unnest(est_tree)
  nodes <- unique(tree$node)
  
  # names(X) <-  paste0("X", 1:ncol(X))
  # data <- data.frame(y, X, group) 
  
  model_mu <- tree %>% 
    dplyr::group_by(node, group) %>% 
    dplyr::summarise(sampled_mu = unique(sampled_mu),
                     sampled_mu_j = unique(sampled_mu_j))
  
  newdata$node <- "root"
  pred <- newdata
  
  
  if(length(nodes) == 1 && nodes == "root"){
    pred_final <-   
      pred %>% dplyr::left_join(
        model_mu, 
        by = c("node", "group"))
    
    return(pred_final)
    
  } else {

  
  results <- dplyr::filter_all(results, any_vars(!is.na(.)))
  results <- results %>% mutate(id = 1:n())
  
  # for each grow followed by prune, remove the parent 
  # to_prune <- dplyr::filter(results, action == 'prune')
  # 
  # 
  # if(nrow(to_prune) > 0){
  #   for(i in 1:nrow(to_prune)){
  #     
  #     # results_change <- results %>% 
  #     #   dplyr::group_by(node) %>% 
  #     #   dplyr::mutate(n_row = 1:n()) %>% 
  #     #   dplyr::ungroup()
  #     
  #     results_change <- results %>% 
  #       filter(id <= to_prune$id[i])
  #     
  #     results_rest <- results %>% 
  #       filter(id > to_prune$id[i])
  #     
  #     results_f <- results_change %>% 
  #       filter(node %in% to_prune$node[i]) %>% 
  #       # getting the closest node
  #       mutate(diff_id = to_prune$id[i] - id) %>% 
  #       dplyr::mutate(to_filter = diff_id %in% c(0, sort(diff_id)[2])) %>% 
  #       dplyr::filter(to_filter) 
  #     
  #     also_to_filter <- results_change %>% 
  #       rowwise() %>% 
  #       mutate(to_filter = str_detect(node, pattern = to_prune$node[i])) %>% 
  #       dplyr::filter(to_filter) 
  #     
  #     results_change <- results_change %>% 
  #       filter(!id %in% c(results_f$id, also_to_filter$id))
  #     
  #     results <- bind_rows(results_change, results_rest)
  #       
  #   } 
  #   
  #   #results <- dplyr::select(results, -to_filter, -n_row)
  # }
  
  # create nodes in new data then left join with main trees
  for(i in 1:nrow(results)){
    pred <- pred %>% 
      dplyr::mutate(
        node = 
          ifelse(
            node == results$node[i], 
            ifelse(!!rlang::sym(results$var[i]) > results$rule[i], 
                   paste(node, results$var[i], "left"), 
                   paste(node, results$var[i], "right")), node))
  }
  # table(tree$node)
  # table(pred$node)
  # results
  
  # antis <- pred %>% dplyr::anti_join(
  #   model_mu, 
  #   by = c("node", "group")) %>% 
  #   dplyr::distinct(node)

  pred_final <-   
    pred %>% dplyr::left_join(
      model_mu, 
      by = c("node", "group"))
  
  return(pred_final)
  
  }
  

} 

#' @name predict_hbm
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Predictions for the B-CART model. 
#' @description This function predicts for the final tree of a 
#' Bayesian CART model. 
#' @param model The model object.  
#' @param newdata The new data to predict. 
#' @return A dataframe with the "prediction" column. 

predict_hbm <- function(model, newdata, formula, group_variable){
  
  model  <- model$final_trees
  response_name <- all.vars(formula)[1]
  group <- dplyr::pull(newdata, !!group_variable)
  m <- stats::model.frame(formula, data = newdata)
  X <- stats::model.matrix(formula, m) %>% 
    as.data.frame()
  
  names(X) <-  paste0("X", 1:ncol(X))
  y <- newdata[ , response_name]
  group <- dplyr::pull(newdata, !!group_variable)
  
  newdata <- data.frame(y, X, group) 
  
  all_results  <- model$results
  P <- length(all_results)
  
  all_preds <- model %>%  
    tibble::add_column(newdata = list(newdata)) %>% 
    dplyr::mutate(final_tree = map(tree_data, ~tail(.x, 1)), 
           preds = pmap(list(results, final_tree, newdata), calculate_preds))  %>% 
    dplyr::select(tree_index, preds)
  
  # all_preds %>% 
  #   tidyr::unnest(preds) %>% 
  #   dplyr::select(tree_index, y, sampled_mu_j, group)  %>% 
  #   View()
  
  #calculate_preds(model$results[[2]], tail(model$tree_data[[2]], 1), newdata)
  
  all_preds_wrangle <- all_preds %>% 
    tidyr::unnest(preds) %>% 
    dplyr::select(tree_index, y, sampled_mu_j, group) %>% 
    #mutate(rn = row_number()) %>% 
    pivot_wider(names_from = tree_index, values_from = sampled_mu_j) %>% 
    unnest() %>% 
    #select(-rn) %>% 
    mutate(pred = rowSums(.[3:ncol(.)])) %>% 
    dplyr::select(y, group, pred)
  

  # for (i in 1:nrow(res)) {
  #   pred <- pred %>% 
  #     dplyr::mutate(node = 
  #                     ifelse(node == 
  #                              res$node[i], ifelse(!!rlang::sym(res$var[i]) > res$rule[i], 
  #                                                  paste(node, res$var[i], "left"), paste(node, res$var[i], 
  #                                                                                         "right")), node))
  # }
  # antis <- pred %>% dplyr::anti_join(model$mu, by = "node") %>% 
  #   dplyr::distinct(node)
  # if (nrow(antis) > 0) {
  #   for (i in 1:nrow(antis)) {
  #     antis$prediction[i] <- model$mu %>%
  #       dplyr::anti_join(pred, by = "node") %>% 
  #       dplyr::filter(stringr::str_detect(node, 
  #                                         antis$node[i])) %>% dplyr::summarise(mu = mean(mu)) %>% 
  #       dplyr::pull(mu)
  #   }
  # }
  # tab_join <- dplyr::bind_rows(model$mu %>% stats::setNames(c("node", 
  #                                                             "prediction")), antis)
  # pred <- pred %>% dplyr::left_join(tab_join, by = "node")
  return(all_preds_wrangle)
}