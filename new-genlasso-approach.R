library(genlasso)
library(CVN)
library(matrixcalc)
library(microbenchmark)
library(glmnet)

# trial for new generalized LASSO estimation
m <- 10 # number of graphs

n1 = .2
n2 = .7
W <- matrix(1, m, m)
D <- create_matrix_D(W, n1, n2)

a = n1^2 + 3*n2^2

A <- diag(a + 1, m) 
A <- diag((n1^2 + 3*n2^2)+0.01, m) 



A - t(D) %*% D
tA = A - t(D) %*% D
eigen(tA)
eigen(diag(n1^2 + 3*n2^2 + 1, m) - t(D) %*% D)

is.positive.semi.definite(A - t(D) %*% D)

y <- rnorm(m, mean = 0, sd = 1)

# apply the generalized LASSO 
out <- genlasso::genlasso(y, diag(1, m), D, minlam = 1)
coef(out, lambda = 1)$beta

#res <- glmnet::glmnet(diag(1,m), y)
#coef(res, s = 30)
microbenchmark::microbenchmark(genlasso::genlasso(y, diag(1, m), D, minlam = 1), 
                               altZ(y, D, diagA = 10, max_iter = 1000), times = 3)


microbenchmark::microbenchmark(altZ(y, D, diagA = n1^2 + 3*n2^2 + 4, max_iter = 1000), times = 3)


altZ(y, D, diagA = n1^2 + 3*n2^2 + 20, max_iter = 1000)


altZ_fast <- function(y, D, A, m, rho = 1, max_iter = 100) { 
  
  p <- nrow(D)
   
  beta_old   <- as.matrix(rep(0, m), ncol = 1)
  beta_new   <- as.matrix(rep(0, m), ncol = 1)
  alpha_new  <- as.matrix(rep(0, max(m, p)), ncol = 1)
  alpha_old1 <- as.matrix(rep(0, max(m, p)), ncol = 1)
  alpha_old2 <- as.matrix(rep(0, max(m, p)), ncol = 1)
  
  A0 <- rho * A 
  A1 <- solve(rho * A + diag(1, m))
  
}


altZ <- function(y, D, diagA = 2, rho = 1, max_iter = 100) { 
  
  m <- length(y)

  A <- diag(diagA, m) 
  
  if (!is.positive.semi.definite(A - t(D) %*% D)) { 
    stop("The matrix A - D'D must be positive semidefinite") 
  }
   
  alpha_new <- as.matrix(rep(0, max(m, nrow(D))), ncol = 1)
  alpha_old1 <- as.matrix(rep(0, max(m, nrow(D))), ncol = 1)
  alpha_old2 <- as.matrix(rep(0, max(m, nrow(D))), ncol = 1)
  
  beta_old <- rep(0, m)
  beta_new <-  rep(0, m)
  #alpha_new <-  rep(0, m)
  #alpha_old1 <-  rep(0, m)
  #alpha_old2 <-  rep(0, m)
  
  #A0 <- rho * A # (1 / rho) * diag(1/3, m)
  #A1 <- solve(rho * A + diag(1, m)) #(1 / rho) * (diag(1/3, m) + diag(1, m))
  

  #tD <- A1 %*% t(D) 
  #A1A0 <- A1 %*% A0 
  #A1y <- A1 %*% y
  #A2 <- A1 %*% t(D)
  c <- 1 / (rho * diagA + 1) 
  c1 <- 1 - c
  cy <- c*y
  cD <- c*D 
  iter <- 0 
  
  repeat{ 
    #beta_new <- A1 %*% (A0 %*% beta_old + y - t(D) %*% (2*alpha_old1 - alpha_old2))
    
    #beta_new <- A1A0 %*% beta_old + A1y - tD %*% (2*alpha_old1 - alpha_old2)
    #beta_new <- beta_old + c*(y - beta_old - tD %*% (2*alpha_old1 - alpha_old2))
    beta_new <- c1*beta_old + cy - crossprod(cD, (2*alpha_old1 - alpha_old2))
    
    #print(iter)
    #print(beta_new)
    
    iter <- iter + 1
    
    #print(sum(abs(beta_new)))
    #print(iter)
    
    alpha_new <- alpha_old1 + rho * D %*% beta_new
    
    alpha_new <- pmax(pmin(alpha_new, 1), -1)
    
    if (sum(abs(beta_new - beta_old)) <= 10^-10 || iter >= max_iter) { 
      break 
    }
    
    beta_old <- beta_new
    alpha_old2 <- alpha_old1
    alpha_old1 <- alpha_new
  }
  
  #print(iter)
  
  beta_new
}

altZ(y, D, max_iter = 1000)
