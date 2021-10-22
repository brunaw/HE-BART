#' @name predict_bcart
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Predictions for the B-CART model. 
#' @description This function predicts for the final tree of a 
#' Bayesian CART model. 
#' @param model The model object.  
#' @param newdata The new data to predict. 
#' @return A dataframe with the "prediction" column. 
predict_bcart <- function(model, newdata){
  
  formula <- model$model_formula
  response_name <- all.vars(formula)[1]
  newdata[, response_name] <- scale(newdata[, response_name]) %>% 
    as.vector()
  m <- stats::model.frame(formula, data = newdata)
  X <- stats::model.matrix(formula, m)
  names(newdata[, colnames(X)]) <- paste("X", 1:length(colnames(X)))
  res <- model$results %>% dplyr::mutate_if(is.factor, as.character)
  newdata$node <- "root"
  pred <- newdata
  for (i in 1:nrow(res)) {
    pred <- pred %>% 
      dplyr::mutate(node = 
                      ifelse(node == 
                               res$node[i], ifelse(!!rlang::sym(res$var[i]) > res$rule[i], 
                                                   paste(node, res$var[i], "left"), paste(node, res$var[i], 
                                                                                          "right")), node))
  }
  antis <- pred %>% dplyr::anti_join(model$mu, by = "node") %>% 
    dplyr::distinct(node)
  if (nrow(antis) > 0) {
    for (i in 1:nrow(antis)) {
      antis$prediction[i] <- model$mu %>%
        dplyr::anti_join(pred, by = "node") %>% 
        dplyr::filter(stringr::str_detect(node, 
                                          antis$node[i])) %>% dplyr::summarise(mu = mean(mu)) %>% 
        dplyr::pull(mu)
    }
  }
  tab_join <- dplyr::bind_rows(model$mu %>% stats::setNames(c("node", 
                                                              "prediction")), antis)
  pred <- pred %>% dplyr::left_join(tab_join, by = "node")
  return(pred)
}