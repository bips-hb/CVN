library(CVN)
data(grid)
m <- 9 # must be 9 for this example

#' Choice of the weight matrix W.
#' (uniform random)
W <- matrix(runif(m*m), ncol = m)
W <- W %*% t(W)
W <- W / max(W)
diag(W) <- 0

# lambdas:
lambda1 = 1:2
lambda2 = 1:2

(cvn <- CVN::CVN(grid, W, lambda1 = lambda1, lambda2 = lambda2, eps = 1e-3, maxiter = 1000, verbose = TRUE))