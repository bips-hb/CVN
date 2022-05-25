
generateX <- function(n = 20, p = 30, mean = 0, sd = 1) { 
  matrix(rnorm(n*p, mean = mean, sd = sd), ncol = p) 
}

X = generateX()

summary(X)
summary(scale(X))
