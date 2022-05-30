library(microbenchmark)

# The sandbox to test ideas and code

# generate a simple dataset X
toy_generateX <- function(n = 20, p = 30, mean = 0, sd = 1) { 
  matrix(rnorm(n*p, mean = mean, sd = sd), ncol = p) 
}

# generate list with square matrices
toy_generateRawDataset <- function(m = 5, p = 30, mean = 0, sd = 1) { 
  lapply(n, function(ni) { 
    toy_generateX(ni, p = p, mean = mean, sd = sd)
  }) 
}

toy_generateRawDataset <- function(m = 5, n = rep(20,m), p = 30, mean = 0, sd = 1) { 
  lapply(n, function(ni) { 
    toy_generateX(ni, p = p, mean = mean, sd = sd)
  }) 
}

# generate a weight matrix 
toy_weightMatrix <- function(m = 5) { 
  matrix(1, nrow = m, ncol = m) 
}



normalized = FALSE

X <- toy_generateRawDataset(mean = 2)

Y_old <- toy_generateRawDataset(mean = 0)
Theta_new <- toy_generateRawDataset()
Theta_old <- toy_generateRawDataset()
Z_new <- toy_generateRawDataset(mean = 0)
  
norm(Y_old[[2]])

# eigen decomposition ------
# m, Theta, Z, Y, Sigma, n_obs, rho = 1, n_cores = 1
m = 10
Sigma <- toy_weightMatrix(m = m)
Theta <- toy_weightMatrix(m = m)
Z <- toy_weightMatrix(m = m)
Y <- toy_weightMatrix(m = m)
rho = 1
n_cores = 1
n_obs = c(100)

A = Sigma - (rho / n_obs[1]) * Z + (rho / n_obs[1]) * Y

eigen_decomposition <- eigen(A)
Q <- eigen_decomposition$vectors 
Lambda <- eigen_decomposition$values

Lambda <- n_obs[1]/(2 * rho) * (-Lambda + sqrt(Lambda^2 + 4*rho / n_obs[1]))
return(Q %*% diag(Lambda) %*% t(Q))

microbenchmark(sum(mcmapply(function(theta_new, theta_old) {norm(theta_new - theta_old)}, 
                            Theta_new, Theta_old)) / sum(mcmapply(norm, Theta_old)), 
               CVN::relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 1), 
               CVN::relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 8))

system.time(sum(mcmapply(function(theta_new, theta_old) {norm(theta_new - theta_old)}, 
       Theta_new, Theta_old)) / sum(mcmapply(norm, Theta_old)))

system.time(
  CVN::relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 1)
)

system.time(
CVN::relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 8)
)



mapply(norm, Theta_new - Theta_old)

mapply('-', Theta_new, Theta_old, SIMPLIFY = FALSE)
Theta_new


lapply()

y = mapply(function(y, theta, z) {y + theta - z}, Y_old, Theta_new, Z_new)

mclapply(1:5, function(i) Y_old[[i]] + Theta_new[[i]] - Z_new[[i]], mc.cores = 8)
b = mapply(function(y, theta, z) {y + theta + z}, Y_old, Theta_new, Z_new, SIMPLIFY = FALSE)

x = lapply(X, function(X) scale(X, center = TRUE, scale = normalized))

Sigma <- lapply(X, cov) 

lapply(X, mean)
lapply(x, mean)

lapply(X, mean)
lapply(X, scale)

r = scale(X[[1]])
mean(scale(X[[1]]))

mapply(function, ...)

X
mu <- lapply(X, mean)

X1 <- mapply('-', X, mu, SIMPLIFY = FALSE)
lapply(X1, mean)
mean(X)
lapply(X, function(x) mean(x))
