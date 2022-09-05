library(genlasso)
library(CVN)
library(matrixcalc)
library(microbenchmark)
library(glmnet)

# trial for new generalized LASSO estimation
m <- 4 # number of graphs

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

altZ(y, D, W, n1, n2, diagA = n1^2 + 3*n2^2 + 20, max_iter = 1000)

microbenchmark::microbenchmark(genlasso::genlasso(y, diag(1, m), D, minlam = 1), 
                               altZ(y, D, W, n1, n2, diagA = n1^2 + 3*n2^2 + 20, max_iter = 1000), times = 3)


#res <- glmnet::glmnet(diag(1,m), y)
#coef(res, s = 30)
microbenchmark::microbenchmark(genlasso::genlasso(y, diag(1, m), D, minlam = 1), 
                               altZ(y, D, diagA = 10, max_iter = 1000), times = 3)


microbenchmark::microbenchmark(altZ(y, D, diagA = n1^2 + 3*n2^2 + 4, max_iter = 1000), times = 3)


altZ(y, D, W, n1, n2, diagA = n1^2 + 3*n2^2 + 20, max_iter = 1000)


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


altZ <- function(y, D, W, eta1, eta2, diagA = 2, rho = 1, max_iter = 100) { 
  
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
    
    x <- crossprod(cD, (2*alpha_old1 - alpha_old2))
    
    alpha <- (2*alpha_old1 - alpha_old2)
    mm <- length(alpha)
    delta <- rep(0, m)
    for (i in 1:m) { 
      #cat(sprintf("i = %d\n", i))
      delta[i] <- eta1 * alpha[i] 
      if (i != m) { 
      for (j in ((i+1):m)) { 
       # cat(sprintf('(%d, %d)\n', i,j))
        #print(W[i,j])
        #print(delta[i])
        #print(alpha[i + j + m - 1])
        
        delta[i] <- delta[i] + eta2*W[i,j]*alpha[i+j+m-1]
      }
      }
      if (i != 1) { 
      for (j in (1:(i-1))) { 
        delta[i] <- delta[i] - eta2*W[j,i]*alpha[i+j+m-1]
      }
      }
    }
    x <- c*delta
    
    beta_new <- c1*beta_old + cy - x
    
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
