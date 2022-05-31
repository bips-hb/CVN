library(CVN) 

load("../MdS/Tests/alladjacs.RData")
load("../MdS/Tests/allsamples.RData")

data <- BigListof_samples[[1]]
true_adj_matrices <- BigListof_adjencies[[1]]

W <- matrix(.5, 9, 9)

res <- CVN(data, W = W, maxiter = 50, lambda1 = 1, lambda2 = 1, rho = 1, verbose = TRUE)

which( abs(res$Theta[[1]]) > 0.1)

which(true_adj_matrices[[1]] == 1)

W

D <- create_matrix_D(W, lambda1 = 1, lambda2 = 1)

y = rep(1,9)
r = genlasso(y, diag(1, 9), D = D, svd = TRUE)

beta = coef(r, lambda = 1)
beta$beta
library(R.utils)
isZero(beta$beta)
beta = beta$beta

beta[which(abs(beta) <= 10^-12)] <- 0

abs(beta$beta) <= 10^-12


