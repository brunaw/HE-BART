#' #' sample_parameters
#' #'
#' #' A function to sample from the main parameter posteriors 
#' #' (tau, mus and mujs)
#' #'
#' #' @name sample_parameters
#' #' @rdname sample_parameters
#' #' @param type The parameter type ("tau", "mu" or "muj")
#' #' @param mu_res_d The mus *** 
#' #' @param N The size of the data
#' #' @param beta A parameter from the tau prior
#' #' @param alpha A parameter from the tau prior
#' #' @param current_tree The current tree
#' #' @param k1 The current value of the k1 parameter
#' #' @param k2 The current value of the k2 parameter
#' #' @param tau_post The current value of the posterior tau
#' #' @param means The current means of the mu_js
#' #' @param M The number of groups
#' #' @param group The group variable
#' #' @param mu_post The current vector of mus
#' #' @export
#' sample_parameters <- function(type = "tau", mu_res_d = NULL,
#'                               N = NULL, beta = NULL, alpha = NULL,
#'                               current_tree = NULL, k1 = NULL,
#'                               tau_post = NULL, 
#'                               k2 = NULL, means = NULL, M = NULL, 
#'                               group = NULL, mu_post = NULL){
#'   
#'   
#'   
#'   if(type == "tau"){
#'     # ones_n <- rep(1, N)
#'     # y <- mu_res_d$y
#'     # Md <- model.matrix(~ factor(mu_res_d$group) - 1)
#'     # MtM <- tcrossprod(x = Md) # Note tcrossprod(X) computes X%*%t(X)
#'     # tMM <- crossprod(x = Md) # crossprod(X) computes t(X)%*%X
#'     # Psi <- k1 * MtM + diag(N)
#'     # W_1 <- Psi + k2 * tcrossprod(x = ones_n)
#'     
#'     errors_y <- mu_res_d %>%
#'       mutate(err = (y - mu_js_sampled)^2) %>%
#'       dplyr::summarise(sum_errors_y = sum(err)) %>%
#'       pull(sum_errors_y)
#'     
#'     alpha_tau <- N/2 + alpha
#'     beta_tau <- errors_y/2 + beta
#'     #beta_tau <-  0.5 * t(y) %*% solve(W_1, y) + beta
#'     sampled <- stats::rgamma(1, alpha_tau, beta_tau)
#'     
#'     
#'   } else if(type == "mu"){
#'     
#'     nam <- unique(current_tree$node)
#'     mu_post <- means <-  sdd <-  vector(length = length(nam))
#'     
#'     for(i in 1:length(nam)){
#'       data_set <- current_tree %>% 
#'         filter(node == nam[i])
#'       M <- model.matrix(~ factor(data_set$group) - 1)
#'       
#'       y <- data_set$y
#'       MtM <- tcrossprod(x = M) # Note tcrossprod(X) computes X%*%t(X)
#'       tMM <- crossprod(x = M) # crossprod(X) computes t(X)%*%X
#'       n <- nrow(data_set)
#'       ones_n <- rep(1, n)
#'       #tau_post <- sqrt(1/tau_post)
#'       
#'       Psi <- k1 * MtM + diag(n)
#'       W_1 <- Psi + k2 * tcrossprod(x = ones_n)
#'       p1_mu <- t(ones_n) %*% solve(Psi, y)
#'       p2_mu <- t(ones_n) %*% solve(Psi, ones_n) + 1 / k2
#'       means[i] <- p1_mu / p2_mu
#'       sdd[i] <- sqrt(1 / tau_post * 1 / p2_mu)
#'       mu_post[i] <- rnorm(1,
#'                           mean = p1_mu / p2_mu,
#'                           sd = sqrt(1 / tau_post * 1 / p2_mu)
#'       )
#'     }
#'     
#'     
#'     
#'     # for(m in 1:n_nodes){
#'     #   means_j <- means[[m]]
#'     #   mean_mu <- (1/k1) * means_j*M/(M/k1 + 1/k2)
#'     #   
#'     #   # Sampling from the posterior distribution of mu ---------
#'     #   var_mu <- (tau_post*(M/k1 + 1/k2))^(-1)
#'     #   
#'     #   mu_post[m] <- stats::rnorm(1, mean = mean_mu, sd = sqrt(var_mu))
#'     # }
#'     
#'     sampled <- mu_post
#'     
#'   } else if(type == "muj"){
#'     n_nodes <- n_distinct(current_tree$node)
#'     mu_js_post <- list()
#'     
#'     njs_mean <-  current_tree %>%
#'       dplyr::group_by(node, group) %>%
#'       summarise(n = n(), mean_y_j = mean(y)) %>% 
#'       arrange(node) %>% 
#'       split(.$node)
#'     
#'     for(m in 1:n_nodes){
#'       nodes_unique <- names(njs_mean)
#'       mean_mu <- c()
#'       var_mu <- c()
#'       
#'       nj_node <- dplyr::pull(njs_mean[[nodes_unique[m]]], n)
#'       y_bar_node <-  dplyr::pull(njs_mean[[nodes_unique[m]]], mean_y_j)
#'       
#'       for(j in sort(unique(group))){
#'         y_bar_j <- y_bar_node[j]
#'         # mean_mu[j] <- ((mu_post[m]/k1) +  y_bar_j * nj_node[j])/(nj_node[j] + 1/k1)
#'         # var_mu[j] <- (tau_post*(nj_node[j] + 1/k1))^(-1)
#'         mean_mu[j] <- (y_bar_j * nj_node[j])/(nj_node[j] + 1/k1)
#'         var_mu[j] <- (2*tau_post*(nj_node[j] + 1/k1))^(-1)
#'       }
#'       mu_js_post[[m]] <- stats::rnorm(M, mean = mean_mu, sd = sqrt(var_mu))
#'     }
#'     sampled <- mu_js_post
#'   }
#'   
#'   return(sampled)
#' }
#' 
#' 
#' # MH_update <- function(
#   current_tree = current_tree, 
#   prior = prior_k1, min_u = min_u, max_u = max_u, pars = pars) {
#   
#   new_k1 <- runif(1, min = min_u, max = max_u) # Or rnorm or whatever else you want to propose
#   current <- calc_all(current_tree, pars$k1, pars)
#   candidate <- calc_all(current_tree, new_k1, pars)
#   
#   if(prior == TRUE){
#     prior_current <- dweibull(pars$k1, shape = 7, 10, log = TRUE)
#     prior_candidate <- dweibull(new_k1, shape = 7, 10, log = TRUE)
#     log.alpha <- candidate - current + prior_candidate - prior_current # ll is the log likelihood, lprior the log prior
#   } else {
#     log.alpha <- candidate - current  
#   }  
#   
#   accept <- log.alpha >= 0 || log.alpha >= log(runif(1))
#   theta <- ifelse(accept, new_k1, pars$k1)
#   return(theta)
# }




#   nam <- unique(current_tree$node)
#   tot_lik <- 0
#   for(i in 1:length(nam)){
#     data_set <- current_tree %>% 
#       filter(node == nam[i])
#     M <- model.matrix(~ factor(data_set$group) - 1)
#     marg <- marg2(y = data_set$y, 
#                        k_1 = new_k1, 
#                        k_2 = pars$k2, 
#                        M = M, 
#                        mu_mu = pars$mu_mu, 
#                        alpha = pars$alpha, 
#                        beta = pars$beta)
#     tot_lik <- marg + tot_lik
#   }
#   
#   return(tot_lik)
#   
# }


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
#
# lk_ratio <- exp_part + first_term







