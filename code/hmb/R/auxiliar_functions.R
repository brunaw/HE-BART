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
  selected_rule <- data %>%
    dplyr::filter(node == sel_node) %>%
    dplyr::mutate(var = !!rlang::sym(variable_index)) %>%
    dplyr::distinct(var) %>%
    dplyr::filter(var > stats::quantile(var, 0.15),
                  var < stats::quantile(var, 0.85)) %>%
    dplyr::pull(var) %>%
    # selecting the cut point
    base::sample(size = 1)

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

  # Extracting the model structure
  mod_str <- stats::model.frame(formula, data = data)

  X <- stats::model.matrix(formula, mod_str)
  # Scaling the response variable
  if(scale_fc == TRUE){
  data[ , response_name] <- scale(data[, response_name]) %>% as.vector()
  }
  group <- data %>% pull(!!group_variable)
  # Defining the response
  y <- data[ , response_name]
  # renaming the covariates
  names(data[ , colnames(X)]) <- paste("X", 1:length(colnames(X)))
  #---------------------------------------------------------------------
  # Initializing accessory columns in the data
  #---------------------------------------------------------------------
  data$node <- "root"           # To save the current node
  data$parent <- "root"         # To save the current parent of each node
  data$d = 0                    # To save the current depth of each node
  data$node_index = 1           # To save the current index of each node
  data$criteria = 'left'        # To initialize the root as a 'left' node
  #---------------------------------------------------------------------

  results <- list(data = data, group = group, y = y, X = X)
  return(results)
}



#' sample_parameters
#'
#' A function to sample from the main parameter posteriors 
#' (tau, mus and mujs)
#'
#' @name sample_parameters
#' @rdname sample_parameters
#' @param type The parameter type ("tau", "mu" or "muj")
#' @param mu_res_d The mus *** 
#' @param N The size of the data
#' @param beta A parameter from the tau prior
#' @param alpha A parameter from the tau prior
#' @param current_tree The current tree
#' @param k1 The current value of the k1 parameter
#' @param k2 The current value of the k2 parameter
#' @param tau_post The current value of the posterior tau
#' @param means The current means of the mu_js
#' @param M The number of groups
#' @param group The group variable
#' @param mu_post The current vector of mus
#' @export
sample_parameters <- function(type = "tau", mu_res_d = NULL,
                              N = NULL, beta = NULL, alpha = NULL,
                              current_tree = NULL, k1 = NULL,
                              tau_post = NULL, 
                              k2 = NULL, means = NULL, M = NULL, 
                              group = NULL, mu_post = NULL){
  

  
  if(type == "tau"){
    # ones_n <- rep(1, N)
    # y <- mu_res_d$y
    # Md <- model.matrix(~ factor(mu_res_d$group) - 1)
    # MtM <- tcrossprod(x = Md) # Note tcrossprod(X) computes X%*%t(X)
    # tMM <- crossprod(x = Md) # crossprod(X) computes t(X)%*%X
    # Psi <- k1 * MtM + diag(N)
    # W_1 <- Psi + k2 * tcrossprod(x = ones_n)

    errors_y <- mu_res_d %>%
      mutate(err = (y - mu_js_sampled)^2) %>%
      dplyr::summarise(sum_errors_y = sum(err)) %>%
      pull(sum_errors_y)
  
    alpha_tau <- N/2 + alpha
    beta_tau <- errors_y/2 + beta
    #beta_tau <-  0.5 * t(y) %*% solve(W_1, y) + beta
    sampled <- stats::rgamma(1, alpha_tau, beta_tau)
    
    
  } else if(type == "mu"){
    
    nam <- unique(current_tree$node)
    mu_post <- means <-  sdd <-  vector(length = length(nam))
    
    for(i in 1:length(nam)){
      data_set <- current_tree %>% 
        filter(node == nam[i])
      M <- model.matrix(~ factor(data_set$group) - 1)
      
      y <- data_set$y
      MtM <- tcrossprod(x = M) # Note tcrossprod(X) computes X%*%t(X)
      tMM <- crossprod(x = M) # crossprod(X) computes t(X)%*%X
      n <- nrow(data_set)
      ones_n <- rep(1, n)
      #tau_post <- sqrt(1/tau_post)
      
      Psi <- k1 * MtM + diag(n)
      W_1 <- Psi + k2 * tcrossprod(x = ones_n)
      p1_mu <- t(ones_n) %*% solve(Psi, y)
      p2_mu <- t(ones_n) %*% solve(Psi, ones_n) + 1 / k2
      means[i] <- p1_mu / p2_mu
      sdd[i] <- sqrt(1 / tau_post * 1 / p2_mu)
      mu_post[i] <- rnorm(1,
                          mean = p1_mu / p2_mu,
                          sd = sqrt(1 / tau_post * 1 / p2_mu)
      )
    }
    
    
    
    # for(m in 1:n_nodes){
    #   means_j <- means[[m]]
    #   mean_mu <- (1/k1) * means_j*M/(M/k1 + 1/k2)
    #   
    #   # Sampling from the posterior distribution of mu ---------
    #   var_mu <- (tau_post*(M/k1 + 1/k2))^(-1)
    #   
    #   mu_post[m] <- stats::rnorm(1, mean = mean_mu, sd = sqrt(var_mu))
    # }
    
    sampled <- mu_post
    
  } else if(type == "muj"){
    n_nodes <- n_distinct(current_tree$node)
    mu_js_post <- list()
    
    njs_mean <-  current_tree %>%
      dplyr::group_by(node, group) %>%
      summarise(n = n(), mean_y_j = mean(y)) %>% 
      arrange(node) %>% 
      split(.$node)

    for(m in 1:n_nodes){
      nodes_unique <- names(njs_mean)
      mean_mu <- c()
      var_mu <- c()
      
      nj_node <- dplyr::pull(njs_mean[[nodes_unique[m]]], n)
      y_bar_node <-  dplyr::pull(njs_mean[[nodes_unique[m]]], mean_y_j)
      
      for(j in sort(unique(group))){
        y_bar_j <- y_bar_node[j]
        # mean_mu[j] <- ((mu_post[m]/k1) +  y_bar_j * nj_node[j])/(nj_node[j] + 1/k1)
        # var_mu[j] <- (tau_post*(nj_node[j] + 1/k1))^(-1)
        mean_mu[j] <- (y_bar_j * nj_node[j])/(nj_node[j] + 1/k1)
        var_mu[j] <- (2*tau_post*(nj_node[j] + 1/k1))^(-1)
      }
      mu_js_post[[m]] <- stats::rnorm(M, mean = mean_mu, sd = sqrt(var_mu))
    }
    sampled <- mu_js_post
  }
  
  return(sampled)
}
