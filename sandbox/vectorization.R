library(CVN)

library(Matrix)
m = 4
p = 10 

set.seed(1)

# generate a list of random matrices
list_of_matrices <- function(m, p) {
  lapply(1:m, function(i) CVN::create_weight_matrix(type = "uniform-random", m = p)) 
}

# turn a list of matrices into a diagonal block matrix
diagonalize <- function(list_matrices) {
  Matrix::bdiag(list_matrices)
}

#Z <- sample_diag_block_matrix(m, p)
#Theta <- sample_diag_block_matrix(m, p)
#Y <- sample_diag_block_matrix(m, p)

# Weight matrix ----------------------------------------------------------------
W <- CVN::create_weight_matrix(type = "uniform-random", m = m)

# Generate data ----------------------------------------------------------------
Theta <- list_of_matrices(m, p)
Y <- list_of_matrices(m, p)

# Old method -------------------------------------------------------------------

Z_old <- updateZ_wrapper(m, p, as.matrix(Theta), as.matrix(Y), W, 
                eta1 = 1, eta2 = 1, a = 100, 
                rho_genlasso = 1, maxiter_genlasso = 100, eps_genlasso = 1e-5, truncate_genlasso = 1e-5)

Z_old <- as.list(Z_old[,1])

# Alternative approach ---------------------------------------------------------
