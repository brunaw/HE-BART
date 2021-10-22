#' @name bcart
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

bcart <- function(formula, dataset, iter = 5000, group_variable, pars,
                  scale = FALSE, p_grow = 0.5, k1_sample = FALSE,
                  min_u = 0, max_u = 20, prior_k1 = FALSE){
  options(dplyr.summarise.inform = FALSE)
  # B-CART Function
  #---------------------------------------------------------------------
  # Handling data
  #---------------------------------------------------------------------
  results_data <- data_handler(formula, dataset, group_variable,
                          scale_fc = scale)
  data <- results_data$data
  names(data)[names(data) == group_variable] <- "group"
  N <- n <- nrow(data)
  group <- results_data$group
  mf <- model.frame(formula, dataset)
  y <- model.extract(mf, "response")
  name_y <- names(mf)[1]
  names(data)[names(data) == name_y] <- "y"
  X <- results_data$X
  
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

  r <- r_k <- vector()          # To save the ratios of grow or prune
  u <- u_k <-  vector()         # To save the sampled uniform values
  sampled_k1 <- vector()
  sampled_k1[1] <- pars$k1
  
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
 
    dn <- max(unique(my_trees[[i]]$d))
    
    alpha_grow <- beta_grow <- 0.5
    p_grow <-   alpha_grow*(1 + dn)^(-beta_grow)
    # Sampling a number from the uniform distribution
    to_do[i] <- stats::runif(1)

    # Checking if the depths of the nodes are bigger than 0
    #     (if 0, we can't prune the tree)
    depths <- sum(unique(my_trees[[i]]$d) > 0)

    # Growing, pruning or staying in the same tree ---------------------
    # Should we grow the tree?
    #p_grow = 0.85
    p_grow = 0.85
    
    if(to_do[i] <  p_grow){
      
      action_taken[i] = "grow"
      # # Selecting the node to grow, uniformly
      drawn_node[i] <- sample(unique(my_trees[[i]]$node), size = 1)
      # # # Selecting the variable and splitting rule, uniformly
      selec_var[i] <- colnames(X)[sample(1:ncol(X), size = 1)]
      # #rule[i]  <- p_rule(variable_index = selec_var[i],
      # #                   data = my_trees[[i]], 
      # #                  sel_node = drawn_node[i])
      # 
      rule[i] <- runif(1, 0.01, 1)
      
      
      # if(change_to_true && i >= 200){
      # 
      #   drawn_node[i] <- "root"
      #   selec_var[i] <- "X1"
      #   rule[i] <- 0.5
      # } 
      # 
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
  
      # Checking whether all of the nodes have the minimum
      # of observations required
      temps <-  my_trees[[i + 1]]  %>%
        dplyr::count(node) %>%
        dplyr::pull(n)
      
      temps_g <-  my_trees[[i + 1]]  %>%
        group_by(node) %>% 
        summarise(n_group = n_distinct(group), 
                  n_group_b = n_group <= 1 , 
                  sum_b = sum(n_group_b))
      

      if(sum(temps <= keep_node) > 0 | sum(temps_g$sum_b) > 0){
        # If all nodes don't have the minimum of observations,
        # just keep the previous tree
        my_trees[[i + 1]] <- my_trees[[i]]
        results[[i+1]] <- results[[i]]
      } else {


        # Saving the parent of the new node to use in the
        # calculation of the transition ratio
          parent_action[i] <- my_trees[[i + 1]] %>%
            dplyr::filter(parent == drawn_node[i]) %>%
            dplyr::distinct(parent) %>%
            dplyr::pull(parent)
        
        results[[i+1]] <- suppressWarnings(
          dplyr::bind_rows(results[[i]],
                           data.frame(
                             node = parent_action[i],
                             var = selec_var[i],
                             rule = rule[i]))
        )

        # Calculating the acceptance ratio for the grow,
        # which uses the ratios of:
        # 1. The transition probabilities of the trees
        # 2. The likelihoods of the trees
        # 3. The probability of the tree structures
        
        r[i] <- ratio_grow(tree = my_trees[[i + 1]],
                           current_node =  parent_action[i],
                           pars = pars, 
                           p = p,
                           current_selec_var = selec_var[i],
                           results_f = results[[i]], 
                           p_grow = p_grow)

      }
      # Should we prune the tree?
    } else if(to_do[i] >  p_grow && depths > 0){
      action_taken[i] = "prune"

      # Selecting node to prune, uniformly
      drawn_node[i] <- sample(unique(my_trees[[i]]$node), size = 1)
      # Detect the two nodes (left and right) to be pruned
      nodes_to_prune <- stringr::str_remove(drawn_node[i], '( right| left)$')
      # Detect the variable that was split in the node that will
      # be now pruned
      variable_in_question <- stringr::str_extract(drawn_node[i],
                                                   'X[0-9][^X[0-9]]*$') %>%
        stringr::str_remove(" left| right")
      # Updating the node that will suffer the prune (returning
      # to the point where it was before the growing)
      new_node <- nodes_to_prune %>% stringr::str_remove('( X[0-9])$')

      my_trees[[i+1]] <- my_trees[[i]] %>%
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


      # Saving the node that was pruned to use in the
      # calculation of the transition ratio
      parent_action[i] <- my_trees[[i]] %>%
        dplyr::filter(node == drawn_node[i]) %>%
        dplyr::distinct(parent) %>%
        dplyr::pull(parent)

      parent_prune <- stringr::str_remove(parent_action[i], '( right| left)$')

      results[[i+1]] <- results[[i]] %>%
        dplyr::filter(!stringr::str_detect(node, parent_prune))

      # Calculating the acceptance ratio for the prune,
      # which uses the ratios of:
      # 1. The transition probabilities of the trees
      # 2. The likelihoods of the trees
      # 3. The probability of the tree structures

      r[i] <- ratio_prune(old_tree = my_trees[[i]],
                          tree = my_trees[[i + 1]],
                          pars = pars,
                          current_node = parent_action[i],
                          var_in_prune = variable_in_question,
                          p = p,
                          nodes_to_prune = nodes_to_prune,
                          results_f = results[[i]], 
                          p_grow = p_grow)

      # Should we stay in the same tree?
    } else {
      action_taken[i] = "same"
      my_trees[[i + 1]] <- my_trees[[i]]
      results[[i+1]] <- results[[i]]
    }

    # Checking if the tree will be accepted or not ---------------------
    if(!identical(my_trees[[i + 1]], my_trees[[i]])){

      # Should we accept the new tree?
      u[i] <- stats::runif(1)

      # Checking if the tree should be accepted or not, based
      # on the acceptance ratio calculated and a value sampled
      # from a uniform distribution
      if(u[i] >= r[i]){
        my_trees[[i + 1]] <- my_trees[[i]]
        results[[i+1]] <- results[[i]]
      }
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

    # ----------------------------------------------------------
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

    #curve(dgamma(x, alpha_tau, beta_tau))

    #tau_post[i + 1] <- 0.4860099

#         errors_y <- suppressWarnings(
#           mu_res[[i]] %>%
#             dplyr::group_by(node) %>%
#             dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
#             dplyr::summarise(sum_errors_y = sum(err_y)))
#
#     err_mus <- suppressWarnings(
#       mu_res[[i]] %>%
#         group_by(node, mu_sampled, mu_js_sampled) %>%
#         slice(1) %>%
#         group_by(node) %>%
#         dplyr::mutate(err_mu = (mu_sampled - mu_js_sampled)^2) %>%
#         dplyr::summarise(sum_errors_mu = sum(err_mu),
#                          mu_sampled = unique(mu_sampled)))
#
#     alpha_tau <- (N + M + 1)/2 + alpha
#     beta_tau <- (mean(errors_y$sum_errors_y)/2 +
#                    beta +
#                    mean(err_mus$sum_errors_mu)/(2 * k1) +
#                    (mean(mu[[i + 1]])^2)/(2 * k2))
#     ### Note: check last term of beta tau !!!!!!!
#     # Is it the mean/sum of mu? How to summarise mu per nodes?
#     tau_post[i + 1] <- stats::rgamma(1, alpha_tau, beta_tau)

    # The version for each node would be:
#     errors_y <- mu_res[[i]] %>%
#       dplyr::group_by(node) %>%
#       dplyr::mutate(err_y = (y - mu_js_sampled)^2) %>%
#       dplyr::summarise(sum_errors_y = sum(err_y))
#
#     err_mus <- mu_res[[i]] %>%
#       group_by(node, mu_sampled, mu_js_sampled) %>%
#       slice(1) %>%
#       group_by(node) %>%
#       dplyr::mutate(err_mu = (mu_sampled - mu_js_sampled)^2) %>%
#       dplyr::summarise(sum_errors_mu = sum(err_mu),
#                        mu_sampled = unique(mu_sampled))
    #   disp_terms <- errors_y %>%
    #     right_join(err_mus, by = "node") %>%
    #     mutate(term = (sum_errors_y)/2 + beta +
    #              (sum_errors_mu)/(2 + k1) + (mu_sampled^2)/(2*k2))
    #
    #   alpha_tau <- (N + M + 1)/2 + alpha
    #   beta_tau <- disp_terms$term
    #
    #   tau_post_all <- stats::rgamma(n_nodes, alpha_tau, beta_tau)

    # New posterior parameters, check if they're correct
    errors_y <- mu_res[[i]] %>%
      mutate(err = (y - mu_js_sampled)^2) %>% 
      dplyr::summarise(sum_errors_y = sum(err)) %>% 
      pull(sum_errors_y)
    
    alpha_tau <- N/2 + alpha
    beta_tau <- errors_y/2 + beta
    tau_post[i + 1] <- stats::rgamma(1, alpha_tau, beta_tau)

    # ------------------------------------------------------------
    if(k1_sample){
    # Sampling K1 -----------------------------------------------
    # We can set these parameters more smartly 
    # new_k1 <- round(runif(n = 1, min = min_u, max = max_u), 2)
    # current_tree = my_trees[[i+1]]
    # 
    # current <- calc_all(current_tree, pars$k1, pars)
    # candidate <- calc_all(current_tree, new_k1, pars)
    # 
    # exp_term <- - current + candidate
    # 
    # if(prior_k1){
    #   # A prior that has high probability in values bigger than 0 
    # exp_term <- - exp_term + (dweibull(new_k1, shape = 7, 10, log = TRUE))
    # }
    # r_k[i] <- min(1, exp(exp_term))
    
    # Should we accept the new tree?
    #u_k[i] <- stats::runif(1)
    
    # Checking if the tree should be accepted or not, based
    # on the acceptance ratio calculated and a value sampled
    # from a uniform distribution
    # if(u_k[i] <= r_k[i]){
    # if(0.6 <= r_k[i]){ 
    #   pars$k1 <- new_k1
    #   sampled_k1[i] <- new_k1
    # } else { 
    #   sampled_k1[i] <- pars$k1
    # }
      
      current_tree = my_trees[[i+1]]
      
      samp_k1 <- MH_update(
        current_tree = current_tree, prior = prior_k1, 
        min_u = min_u, max_u = max_u, pars = pars)
      
      if(samp_k1 == pars$k1){
        sampled_k1[i] <- pars$k1
      } else {
        pars$k1 <- samp_k1
        sampled_k1[i] <- samp_k1
      }
      
   
    }
    
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
    results = results[[i]],
    #mu = mu_res[[i-1]],
    mus = mu,
    model_formula = formula,
    selec_var = selec_var,
    rule = rule,
    sampled_k1 = sampled_k1, 
    u_k = u_k, 
    r_k = r_k
  ))
}

