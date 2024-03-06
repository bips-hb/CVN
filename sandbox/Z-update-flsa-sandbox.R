library(wflsa)

m = 100
p = 10

W <- CVN::create_weight_matrix('full', m = m)

connList <- CVN::create_connListObj(W)


Theta <- lapply(1:m, function(i) {
  mat <- CVN::create_weight_matrix('uniform-random', m = p)
  diag(mat) <- rnorm(p)
  mat
})

Y <- lapply(1:m, function(i) {
  mat <- CVN::create_weight_matrix('uniform-random', m = p)
  diag(mat) <- rnorm(p)
  mat
})

# update the diagonal
Z <- lapply(1:m, function(i) {
  diag(diag(Theta[[i]]) - diag(Y[[i]]))}
)
 
# update the off-diagonal
co <- combn(p,2)

apply(co, 2, function(pair) {
  s <- pair[1]
  t <- pair[2]
  print(c(s,t))
  
  y <- sapply(1:m, function(i) Theta[[i]][s,t] + Y[[i]][s,t])
  if(connList$connList_is_null) {
    fit <- flsa::flsa(y, lambda1 = eta1, lambda2 = 0)
  } else {
    fit <- flsa::flsa(y, lambda1 = eta1, lambda2 = eta2, connListObj = connList$connList)
  }
  
  fit <- as.vector(fit)
  sapply(1:m, function(i) {
    Z[[i]][s,t] <<- fit[i]
    Z[[i]][t,s] <<- fit[i]
  })
})

s <- 2
t <- 5

eta1 = .001
eta2 = 1

y <- sapply(1:m, function(i) Theta[[i]][s,t] + Y[[i]][s,t])


if(connList$connList_is_null) {
  fit <- as.vector(flsa::flsa(y, lambda1 = eta1, lambda2 = 0))
} else {
  fit <- as.vector(flsa::flsa(y, lambda1 = eta1, lambda2 = eta2, connListObj = connList$connList))
}

as.vector(wflsa::wflsa(y, W, lambda1 = eta1, lambda2 = eta2))

microbenchmark::microbenchmark(
  as.vector(flsa::flsa(y, lambda1 = eta1, connListObj = connList$connList)),
  as.vector(wflsa::wflsa(y, W = CVN::create_weight_matrix('uniform-random', m), lambda1 = eta1, lambda2 = eta2))
)

as.vector(fit)




