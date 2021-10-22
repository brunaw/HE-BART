#' @name bcart_fixed
#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}.
#' @export
#' @title B-CART model.
#' @description This function runs a BCART model and returns the tree and
#' other results obtained in the last iteration of the MCMC.
#' @param formula The model formula.
#' @param data The data to be used in the modeling.
#' @param iter The number of iterations for the MCMC.
#' @param group_variable The grouping variable.
#' @param pars The hyperparameter set
#' @param scale Logical to decide whether to scale y or not
#' @return A list containing:
#'  sigma^2_y, the current errors, the final tree, the
#' ratios used in the MCMC step and the uniform values sampled,
#  the action taken at each step, the "path" go the final tree,
#' the final vector of mus and the model formula.
#' @details
#' Priors used ----------------------------------------------------------
#' sigma^2 ~ InvGamma(nu/2, lambda*nu/2), with
#' the parameters chosen in a way that there is a high probability mass
#' 'before' the RMSE of a linear model
#'
#' mu_i ~ Normal(mu_mu, sd_mu^2), with
#' mu_mu = 0
#' sd_<- = (max(y))/2
#'
#' Posteriors (both conjugate) ------------------------------------------
#' sigma^2 ~ InvGamma((nu + n) / 2, 1/2 *(sum_i error^2_i + lambda * nu)))
#' Note: the error is calculated with the current mu_i for each node
#'
#' mu_i ~ Normal(mu_post, sigma_mu_post), with
#' sigma_mu_post = 1/(1/sd_mu + n_in_the_node/sigma^2)
#' mu_post = (0 + (n_in_the_node/sigma^2) * node.avgResponse* sigma_mu_post
#' n_in_the_node = observations in each node
#' node.avgResponse = average of the response in each node

bcart_fixed <- function(formula, dataset, iter = 5000, group_variable, pars,
                  scale = FALSE, p_grow = 0.5){
  options(dplyr.summarise.inform = FALSE)
  # B-CART Function
  #---------------------------------------------------------------------
  # Handling data
  #---------------------------------------------------------------------
  results <- data_handler(formula, dataset, group_variable,
                          scale_fc = scale)
  data <- results$data
  N <- n <- nrow(data)
  group <- results$group
  y <- data$y
  X <- results$X
  #---------------------------------------------------------------------
  # Defining current distribution parameters
  #---------------------------------------------------------------------
  # Prior hyperparameters -------------
  M <- n_distinct(group)
  beta <- pars$beta
  alpha <- pars$alpha
  mu_mu <- pars$mu_mu
  k1 <- pars$k1
  k2 <- pars$k2
  
  # Priors for tau, tau_mu and tau_mu_j  -------------
  #tau_pr <- stats::rgamma(n = 1, 1/alpha, beta)
  tau_pr <- stats::rgamma(n = 1, 1/alpha, beta)
  tau_mu_j <-  k1 * (1/tau_pr)
  tau_mu <-  k2 * (1/tau_pr)
  
  #alpha_tau <- N/2 + alpha
  #alpha_tau <- 2/2 + alpha
  #M_mat <- stats::model.matrix(y ~ factor(group) - 1, data)
  #psi <- (k1 * M_mat %*% t(M_mat)) + diag(N)
  
  #W_0 <- rep(mu_mu, N)
  #W_1 <- (k2 * diag(x = 1, nrow = N, ncol = 1) %*%
  # t(diag(x = 1, nrow = N, ncol = 1))) + psi
  #in_W_1 <- solve(W_1)
  #beta_tau <- (t((y - W_0)) %*% in_W_1 %*% (y - W_0))/2 + beta
  
  # Priors for all the mus  -------------
  #mu_main <-  stats::rnorm(n = 1, 0, sd = sqrt(tau_mu))
  mu_main <- mean(y)
  mu_j <- data %>%
    group_by(group) %>%
    summarise(m = mean(y)) %>% pull(m)
  
  # Creating objects to save the results -------------
  # Minimum batch for each node
  keep_node <- 0.05 * nrow(data)
  p <- ncol(X)
  to_do <- vector()         # Number of actions that can be taken
  
  #---------------------------------------------------------------------
  # Initializing useful vectors
  #---------------------------------------------------------------------
  # For grow or prune
  action_taken = vector()       # To save the actions taken in the algorithm
  selec_var = vector()          # To save the selected variable when growing
  rule = vector()               # To save the selected splitting rule when growing
  drawn_node = vector()         # To save the selected node to grow or prune
  
  r <- vector()              # To save the ratios of grow or prune
  u <- vector()              # To save the sampled uniform values
  
  # For the trees  -------------
  my_trees <-  list()            # To save each new tree
  my_trees[[1]] <- data         # Initializing the first tree as the
  # data will be the  'root' tree, with no nodes
  
  # For the sampling of posterior values  -------------
  tau_post <- vector()            # To save posterior values of
  mu <- mu_js <- list()           # tau_mu, mu and mu_js
  mu[[1]] <- mu_main
  mu_js[[1]] <- mu_j
  tau_post[1] <-  tau_pr
  
  errors <- list()          # To save the value of the errors used in
  parent_action <- vector()   # the posterior of sigma
  results <- list()
  results[[1]] <- data.frame(node = NA, var = NA, rule = NA)
  mu_res <- list()
  #---------------------------------------------------------------------
  # A simple progress bar
  pb <- progress::progress_bar$new(
   format = "  Iterations of the BCART model [:bar]  :current/:total (:percent)",
   clear = FALSE, width = 60, total = iter)

  # Adding first mu and mu_j values
  add_mu_js <- my_trees[[1]] %>%
    dplyr::distinct(node, group) %>%
    dplyr::arrange(node, group) %>%
    dplyr::mutate(mu_js_sampled = mu_j)
  
  mu_res[[1]] <- my_trees[[1]] %>%
    dplyr::distinct(node) %>%
    dplyr::arrange(node) %>%
    dplyr::mutate(mu_sampled = mu[[1]]) %>%
    dplyr::right_join(my_trees[[1]], by = "node") %>%
    dplyr::right_join(add_mu_js, by = c("node", "group"))
  
  my_trees[[1]] <- mu_res[[1]]
  
  
  # Loop to perform the BCART model
  for(i in 1:iter){
    
    #dn <- max(unique(my_trees[[i]]$d))
    
    # Checking if the depths of the nodes are bigger than 0
    #     (if 0, we can't prune the tree)
    #depths <- sum(unique(my_trees[[i]]$d) > 0)
    
    # Growing, pruning or staying in the same tree ---------------------
    # Should we grow the tree?

    if(i == 1){
      rule[i]   <- 0.5
      drawn_node[i] <- "root"
      selec_var[i] <- "X1"
    } else if(i == 2){
      rule[i]   <- 0.5
      drawn_node[i] <- "root X1 right"
      selec_var[i] <- "X2"
    }
    
    if(i < 3){
      my_trees[[i+1]] <- my_trees[[i]] %>%
        dplyr::mutate(
          # Increasing the depth of the node giving the grow
          d =  ifelse(node == drawn_node[i], d + 1, d),
          # Updating the parent of the splitted node
          parent = ifelse(node == drawn_node[i], drawn_node[i], parent),
          
          # Changing the node "side" of each observation: left and right
          criteria = ifelse(
            node == drawn_node[i],
            ifelse(!!rlang::sym(selec_var[i]) > rule[i],
                   "left", "right"), "no split"),
          
          # Updating the node accordingly to the new split
          node = ifelse(node == drawn_node[i],
                        ifelse(!!rlang::sym(selec_var[i]) > rule[i],
                               paste(node, selec_var[i], "left"),
                               paste(node, selec_var[i], "right")), node),
          # Updating the node index
          node_index =  as.numeric(as.factor(node)))
    } else if(i >= 3){
      rule[i]   <- 0.5
      drawn_node[i] <- "root X1 left"
      selec_var[i] <- "X2"
      
      my_trees[[i+1]] <- my_trees[[3]]
    }
    
    
    # Updating posteriors -----------
    n_nodes <- my_trees[[i + 1]] %>%
      dplyr::distinct(node) %>%
      nrow()
    
    # Creating a vector to save a mu for each node
    mu_post <- vector(length = n_nodes)
    mu_js_post <- list()
    
    # Sampling from the posterior distribution of mu_j -----------
    # For each node, need to get mu from the parent it seems
    njs <-  my_trees[[i + 1]] %>%
      dplyr::group_by(node, group) %>%
      dplyr::count()
    
    split_nodes <- njs %>% split(.$node)
    n_j_means <- my_trees[[i+1]] %>%
      group_by(node, group) %>%
      summarise(mean_y_j = mean(y)) %>%
      arrange(group) %>%
      split(.$node)
    
    for(m in 1:n_nodes){
      nodes_unique <- unique(njs$node)
      mean_mu <- c()
      var_mu <- c()
      
      nj_node <- split_nodes[[nodes_unique[m]]] %>% dplyr::pull(n)
      y_bar_node <-  n_j_means[[nodes_unique[m]]] %>% dplyr::pull(mean_y_j)
      
      for(j in sort(unique(group))){
        y_bar_j <- y_bar_node[j]
        mean_mu[j] <- ((mu_post[m]/k1) +  y_bar_j * nj_node[j])/(nj_node[j] + 1/k1)
        var_mu[j] <- (tau_post[i]*(nj_node[j] + 1/k1))^(-1)
      }
      mu_js_post[[m]] <- stats::rnorm(M, mean = mean_mu, sd = sqrt(var_mu))
    }
    
    mu_js[[i + 1]] <- mu_js_post
  
    
    means <- lapply(mu_js_post, mean)
    
    for(m in 1:n_nodes){
      #means_j <- means$mean_mu_j[m]
      means_j <- means[[m]]
      mean_mu <- (1/k1) * means_j*M/(M/k1 + 1/k2)
      
      # Sampling from the posterior distribution of mu ---------
      var_mu <- (tau_post[i]*(M/k1 + 1/k2))^(-1)
      
      mu_post[m] <- stats::rnorm(1, mean = mean_mu, sd = sqrt(var_mu))
    }
    mu[[i + 1]] <- mu_post
    
    
    # -----------------------------------------------------------
    # Update these values in the tree -------------
    # -----------------------------------------------------------
    if(!(i == 1)){
      add_mu_js <- my_trees[[i+1]] %>%
        dplyr::distinct(node, group) %>%
        dplyr::arrange(node, group) %>%
        dplyr::mutate(mu_js_sampled = unlist(mu_js_post))
      
      mu_res[[i]] <- my_trees[[i+1]] %>%
        dplyr::distinct(node) %>%
        dplyr::arrange(node) %>%
        dplyr::mutate(mu_sampled = mu[[i+1]])
      
      mu_res[[i]] <- add_mu_js %>%
        right_join(mu_res[[i]], by = "node")
      
      mu_res[[i]] <- mu_res[[i]] %>%
        dplyr::left_join(
          select(my_trees[[i+1]], -mu_js_sampled, -mu_sampled),
          by = c("node", "group"))
      
      my_trees[[i + 1]] <- mu_res[[i]]
    }
    
    # -----------------------------------------------------------
    # Sampling from the posterior distribution of tau -----------

    # New posterior parameters, check if they're correct
    errors_y <- mu_res[[i]] %>%
      dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
      dplyr::summarise(sum_errors_y = sum(err_y))
    
    alpha_tau <- N/2 + alpha
    beta_tau <- errors_y$sum_errors_y/2 + beta
    
    #alpha_tau/beta_tau
    #sqrt(1/(alpha_tau/beta_tau))
    
    # Issue: beta is much higher than it should be 
    tau_post[i + 1] <- stats::rgamma(1, alpha_tau, beta_tau)

    # ------------------------------------------------------------
    pb$tick()
  }
  
  cat("
# ---------------------------------------------
The tree has finished growing, with a final final posterior precision:", tau_post[i-1], "\nThe minimum node depth is", min(my_trees[[i-1]]$d), "and the maximum is ",
      max(my_trees[[i-1]]$d),
      '\n# ---------------------------------------------')
  
  return(list(
    tau_post = tau_post,
    #errors = errors,
    final_tree = my_trees[[i-1]],
    ratios = r,
    samp_unif = u,
    nodes = drawn_node,
    action_taken = action_taken,
    trees = my_trees,
    #results = results[[i]],
    #mu = mu_res[[i-1]],
    mus = mu,
    model_formula = formula,
    selec_var = selec_var,
    rule = rule
  ))
}

