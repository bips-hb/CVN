library(CVN)
data(grid)
m <- 9 # must be 9 for this example

#' Choice of the weight matrix W.
#' (uniform random)
W <- matrix(runif(m*m), ncol = m)
W <- W %*% t(W)
W <- W / max(W)
diag(W) <- 0

W <- CVNSim::create_weight_matrix("grid")

# lambdas:
lambda1 = c(.1, 1:4)
lambda2 = c(.1, 1:4)

(cvn <- CVN::CVN(grid, W, n_cores = 7, warmstart = TRUE, lambda1 = lambda1, lambda2 = lambda2, verbose = TRUE))
