library(genlasso)
library(CVN)
library(matrixcalc)
library(microbenchmark)
library(glmnet)

# trial for new generalized LASSO estimation
m <- 8 # number of graphs

global_rho <- 1
lambda1 <- .2
lambda2 <- .7

n1 <- global_rho * lambda1
n2 <- global_rho * lambda2
W <- matrix(1, m, m)
D <- create_matrix_D(W, lambda1, lambda2, global_rho)
c <- nrow(D)

a = n1^2 + 3*n2^2

y <- rnorm(m, mean = 0, sd = 1)

# apply the generalized LASSO 
out <- genlasso::genlasso(y, diag(1, m), D, minlam = 1)
coef(out, lambda = 1)$beta

altZ(y, D, W, n1, n2, diagA = 20, rho = 1, max_iter = 1000)
  

CVN::aug_genlasso(y, W, m, c, lambda1, lambda2, global_rho, a, rho = 1, max_iter = 1000, eps = 10^-10)
  