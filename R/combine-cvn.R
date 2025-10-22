#' Combine CVNs of size 1
#' 
#' Function to combine single CVNs. Handy, if CVNs are computed in parallel using
#' foreach loop and n_cores = 1
#'
#' @param cvn_list a list of CVNs, each of size 1 (one lambda/gamma pair)
#' @param minimal  logical; stripped CVN (default = TRUE)
#'
#' @returns a CVN of size length(cvn_list)
#' @export
#'
#' @examples
#' data(grid)
#'
#' #' Choice of the weight matrix W. Each of 2 covariates has 3 categories
#' #' (uniform random)
#' W_grid <- create_weight_matrix("grid", k = 3, l = 3)
#' 
#' fit1 <- CVN(data = grid, W = W_grid, lambda1 = 1, lambda2 = 0.5, 
#'             gamma_ebic = 0.75,
#'             n_cores = 1, eps = 1e-2, maxiter = 200, minimal = TRUE)
#'             
#' fit2 <- CVN(data = grid, W = W_grid, lambda1 = 2, lambda2 = 0.75, 
#'             gamma_ebic = 0.75,
#'             n_cores = 1, eps = 1e-2, maxiter = 200, minimal = TRUE)
#'             
#' (fit <- combine_cvn(list(fit1, fit2)))            

combine_cvn <- function(cvn_list, minimal = TRUE){

  # cvn_list is a list of cvn fits with one lambda-value combination
  is_cvn <- sapply(cvn_list, is)
  stopifnot("all list elements must be of type 'cvn'" = all(is_cvn == is_cvn[1])) 
  
  N_lambda <- sapply(cvn_list, function(fit) fit$n_lambda_values)
  stopifnot("all CVNs must be minimal" = all(N_lambda == 1))
  
  is_minimal <- sapply(cvn_list, function(fit) fit$minimal)
  stopifnot("all CVNs must be minimal" = all(is_minimal == is_minimal[1])) 
  
  m_values <- sapply(cvn_list, function(fit) fit$m)
  stopifnot("m must be equal in all CVNs" = all(m_values == m_values[1])) 
  
  p_values <- sapply(cvn_list, function(fit) fit$p)
  stopifnot("p must be equal in all CVNs" = all(p_values == p_values[1])) 
  
  maxiter_values <- sapply(cvn_list, function(fit) fit$maxiter)
  stopifnot("maxiter must be equal in all CVNs" = all(maxiter_values == maxiter_values[1])) 
  
  rho_values <- sapply(cvn_list, function(fit) fit$rho)
  stopifnot("rho must be equal in all CVNs" = all(rho_values == rho_values[1])) 
  
  eps_values <- sapply(cvn_list, function(fit) fit$eps)
  stopifnot("eps must be equal in all CVNs" = all(eps_values == eps_values[1]))   
  
  truncate_values <- sapply(cvn_list, function(fit) fit$truncate)
  stopifnot("truncate must be equal in all CVNs" = all(truncate_values == truncate_values[1])) 
  
  rho_genlasso_values <- sapply(cvn_list, function(fit) fit$rho_genlasso)
  stopifnot("rho_genlasso must be equal in all CVNs" = 
              all(rho_genlasso_values == rho_genlasso_values[1])) 
  
  eps_genlasso_values <- sapply(cvn_list, function(fit) fit$eps_genlasso)
  stopifnot("eps_genlasso must be equal in all CVNs" = 
              all(eps_genlasso_values == eps_genlasso_values[1]))  
  
  maxiter_genlasso_values <- sapply(cvn_list, function(fit) fit$maxiter_genlasso)
  stopifnot("maxiter_genlasso must be equal in all CVNs" = 
              all(maxiter_genlasso_values == maxiter_genlasso_values[1])) 
  
  truncate_genlasso_values <- sapply(cvn_list, function(fit) fit$truncate_genlasso)
  stopifnot("truncate_genlasso must be equal in all CVNs" = 
              all(truncate_genlasso_values == truncate_genlasso_values[1]))   
  
  normalized_values <- sapply(cvn_list, function(fit) fit$normalized)
  stopifnot("normalized must be equal in all CVNs" = 
              all(normalized_values == normalized_values[1]))
  
  gamma_values <- sapply(cvn_list, function(fit) fit$gamma_ebic)
  if(any(gamma_values != gamma_values[1])){warning("gamma for ebic NOT equal in all CVNs!")}
  
  W_mat <- sapply(cvn_list, function(fit) fit$W)
  all_equal <- apply(W_mat, 2, function(col) all(col == W_mat[,1]))
  stopifnot("W must be equal in all CVNs" = all(all_equal)) 
  
  N <- sapply(cvn_list, function(fit) fit$n_obs)
  all_equal <- apply(N, 2, function(col) all(col == N[,1]))
  stopifnot("n_obs must be equal in all CVNs" = all(all_equal))   
  
  # combine
  combined_cvn <- cvn_list[[1]]
  combined_cvn$adj_matrices <- lapply(cvn_list, function(x) x$adj_matrices[[1]])
  combined_cvn$warmstart <- sapply(cvn_list, function(fit) fit$warmstart)
  combined_cvn$hits_border_aic <- sapply(cvn_list, function(fit) fit$hits_border_aic)
  combined_cvn$hits_border_bic <- sapply(cvn_list, function(fit) fit$hits_border_bic)
  tmp <- as.data.frame(do.call(rbind, lapply(cvn_list, function(x) x$results)))
  if(any(gamma_values != gamma_values[1])){
    tmp$gamma_ebic <- gamma_values
  }
  combined_cvn$results <- tmp[order(-tmp$lambda1, -tmp$lambda2),]
  combined_cvn$results$id <- 1:nrow(combined_cvn$results)
  combined_cvn$n_lambda_values <- nrow(combined_cvn$results)  

  if ('hits_border_aic' %in% names(combined_cvn)) {
    combined_cvn <- within(combined_cvn, rm(hits_border_aic))
  }
  if ('hits_border_bic' %in% names(combined_cvn)) {
    combined_cvn <- within(combined_cvn, rm(hits_border_bic))
  }
  
  # Set the class back to "irgendwas"
  class(combined_cvn) <- c("cvn", "list")
  return(combined_cvn)
  
}
