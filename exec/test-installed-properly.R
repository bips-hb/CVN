library(CVN)

cat("Running a simple example...\n\n")

cat("Applying CVN to the data...\n")

data(grid)
m <- 9 # must be 9 for this example

#' Choice of the weight matrix W.
#' (uniform random)
W <- matrix(runif(m*m), ncol = m)
W <- W %*% t(W)
W <- W / max(W)
diag(W) <- 0

# lambdas:
lambda1 = c(1, 2)
lambda2 = c(1, 2)

cvn <- CVN::CVN(grid, W, warmstart = TRUE, eps = 1e-4, maxiter = 1000,
                 lambda1 = lambda1, lambda2 = lambda2, verbose = TRUE)

cat("DONE CVN has been estimated\n")
cat("Trying plotting functions...\n")
plot_hamming <- CVN::plot_hamming_distances_cvn(cvn)
plot_aic <- CVN::plot_aic(cvn)
plot_cvn <- CVN::visnetwork_cvn(cvn)

plot_cvn$plots[[1]][[1]]
cat("DONE creating plots\n")

