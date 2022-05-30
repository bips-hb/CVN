# The sandbox to test ideas and code

# generate a simple dataset X
toy_generateX <- function(n = 20, p = 30, mean = 0, sd = 1) { 
  matrix(rnorm(n*p, mean = mean, sd = sd), ncol = p) 
}

# generate raw dataset
toy_generateRawDataset <- function(m = 5, n = rep(20,m), p = 30, mean = 0, sd = 1) { 
  lapply(n, function(ni) { 
    toy_generateX(ni, p = p, mean = mean, sd = sd)
  }) 
}

# generate a weight matrix 
toy_weightMatrix <- function(m = 5) { 
  matrix(rep(1, m*m), nrow = m) 
}


normalized = FALSE

X <- toy_generateRawDataset(mean = 2)

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
