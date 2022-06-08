library(CVN) 
library(CVNSim)
library(parallel)
library(profvis)
library(quadprog)
library(LowRankQP)
library(microbenchmark)

set.seed(1)

m = 4
n = 100
p = 5

graphs <- lapply(1:m, function(i) CVNSim::generate_graph(p = p))
data <- lapply(graphs, function(g) generate_data_given_adjaceny_matrix(n = n, g$adj_matrix))

m <- length(data)       # total number of graphs  
p <- ncol(data[[1]])    # total number of variables
n_obs <- sapply(data, function(X) nrow(X))    # no. of observations per graph 


Sigma <- lapply(1:m, function(i) cov(data[[i]])*(n_obs[i] - 1) / n_obs[i]) 

W <- matrix(1, m, m)


est <- huge(data[[1]], method = "glasso", verbose = FALSE)
q = huge.select(est)
  
microbenchmark(CVN::CVN(data, W = W, maxiter = 1000, epsilon = 1e-2, lambda1 = 20, lambda2 = 2, rho = 1, n_cores = 8, warmstart = FALSE), 
               CVN::CVN(data, W = W, maxiter = 1000, epsilon = 1e-2, lambda1 = 20, lambda2 = 2, rho = 1, n_cores = 8), times = 5)

D <- as.matrix(res$D)
y <- rnorm(m)

library(JGL)
j = JGL(data,penalty="fused",.1,.1,rho=1,weights="equal",penalize.diagonal=FALSE,
    maxiter=500,tol=1e-5,warm=NULL,return.whole.theta=FALSE, screening="fast",
    truncate = 1e-5)

cvn = CVN::CVN(data, W = W, maxiter = 1000, verbose = TRUE, epsilon = 1e-5, lambda1 = 3, lambda2 = .1, rho = 1, n_cores = 8, warmstart = FALSE) 



ROwnsOptim <- function(y, D) { 

  fn <- function(u, y, D) {
    .5 * t((t(D) %*% u - y))%*%(t(D) %*% u - y)
  }

  est = optim(par = as.matrix(rep(0 , nrow(D))), fn, 
            method = "L-BFGS-B", 
            lower = rep(-1,nrow(D)), 
            upper = rep(1,nrow(D)), 
            y = y, 
            D = D)
 
  b = y - t(D) %*% est$par
  b[abs(b) <= 1e-10] <- 0
  return(b)
}

est = ROwnsOptim(y, D)

GENLASSO <- function(y, D) { 
  out <- genlasso(y, diag(1, m), D, minlam = 1)
  beta <- coef(out, lambda = 1)$beta
  beta[abs(beta) <= 1e-10] <- 0
  beta
}

GENLASSO(y, D)



OPTIMX <- function(y, D) { 
  
  fn <- function(u, y, D) {
    .5 * t((t(D) %*% u - y))%*%(t(D) %*% u - y)
  }
  
  m = length(y)
  
  # est = optimx::optimr(par = as.matrix(rep(0 , nrow(D))), fn, 
  #                   method = "L-BFGS-B", 
  #                   lower = rep(-1,nrow(D)), 
  #                   upper = rep(1,nrow(D)), 
  #                   y = y, 
  #                   D = D)
  
  est = optimx::opm(par = as.matrix(rep(0 , nrow(D))), fn, 
                       method = "L-BFGS-B", 
                       lower = rep(-1,nrow(D)), 
                       upper = rep(1,nrow(D)), 
                       control = list(
                         trace = 0,
                         save.failures = FALSE, 
                         all.methods = FALSE,
                         kkt = FALSE
                       ), 
                       y = y, 
                       D = D)
  
  #print(as.matrix(est[1:nrow(D)]))
  #print(t(D))
  
  b = y - t(D) %*% t(as.matrix(est[1:nrow(D)]))
  b[abs(b) <= 1e-10] <- 0
  
  #b = y - t(D) %*% est$par
  #b[abs(b) <= 1e-10] <- 0
  return(b)
}

est = OPTIMX(y, D)

microbenchmark(ROwnsOptim(y,D), GENLASSO(y, D), OPTIMX(y,D))
