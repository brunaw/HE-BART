#' @name p_rule
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title Rule selection.
#' @description Selects a value to split the tree in a grow step.
#' @param variable_index The variable to create the split.
#' @param data The current tree.
#' @param sel_node The node to break from.
#' @return The selected splitting value.

p_rule <- function(variable_index, data, sel_node){
  var <- data %>%
    dplyr::filter(node == sel_node) %>%
    dplyr::pull(variable_index) %>% 
    unique()
  
  # selecting the cut point
  # selected_rule <- sample(var[
  #   var > stats::quantile(var, 0.05) & var < stats::quantile(var, 0.95)], 
  #   size = 1)
  
  
  if(length(var) < 2){
    return(NA)
  } else {
    var <- var[!c(var %in% c(min(var), max(var)))]
    selected_rule <- sample(var, size = 1)
    
  }
  # if(selected_rule %in% c(min(var), max(var))){
  #   selected_rule <- sample(var, size = 1)  
  # }
  
  
  return(selected_rule)
}

#' Pipe operator
#'
#' See \code{\link[magrittr]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
NULL

#' Model data handler
#'
#' A function to adjust the data
#'
#' @name data_handler
#' @rdname data_handler
#' @param formula The model formula.
#' @param data The modelling dataset.
#' @param group_variable The grouping variable.
#' @param scale_fc Logical to decide whether to scale y or not
#' @keywords internal
#' @export

data_handler <- function(formula, data, group_variable, scale_fc){
  #---------------------------------------------------------------------
  # Extracting the data and the response from the formula
  # --------------------------------------------------------------------
  # Removing the intercept
  formula <- stats::as.formula(paste(c(formula), "- 1"))
  response_name <- all.vars(formula)[1]
  names_x <- all.vars(formula[[3]])

  data <- dplyr::select(data, c(!!response_name, !!names_x, !!group_variable))
  # Extracting the model structure
  mod_str <- stats::model.frame(formula, data = data)

  X <- stats::model.matrix(formula, mod_str) %>% 
    as.data.frame()
  # Scaling the response variable
  if(scale_fc == TRUE){
  data[ , response_name] <- scale(data[, response_name]) %>% as.vector()
  }
  group <- dplyr::pull(data, !!group_variable)
  # Defining the response
  y <- data[ , response_name]
  # renaming the covariates
  depara_names <- tibble(original = names(X),
                         new = paste0("X", 1:ncol(X)))
  
  
  names(X) <-  paste0("X", 1:ncol(X))
  data <- data.frame(y, X, group) 
  
  #---------------------------------------------------------------------
  # Initializing accessory columns in the data
  #---------------------------------------------------------------------
  data$node <- "root"           # To save the current node
  data$parent <- "root"         # To save the current parent of each node
  data$d = 0                    # To save the current depth of each node
  data$node_index = 1           # To save the current index of each node
  data$criteria = 'left'        # To initialize the root as a 'left' node
  #---------------------------------------------------------------------

  results <- list(data = data, group = group, y = y, X = X, names = depara_names)
  return(results)
}



