library(quadprog)

m = 4
fn <- function(u, y, D = diag(1,m)) { 
  print(y)
   .5 * t((t(D) %*% u - y))%*%(t(D) %*% u - y) 
}

x = optim(par = as.matrix(seq(.5,m)), fn, method = "L-BFGS-B", lower = rep(-1,m), upper = rep(1,m), y = rep(1:m))


D <- res$D
m = res$m

fn <- function(u, y, D) { 
  .5 * t((t(D) %*% u - y))%*%(t(D) %*% u - y) 
}

t(D) %*% u
dim(D)
dim(t(D))
dim(t(u))

r = m*(m-1)/2

est = optim(par = rep(0,r+m), fn, method = "L-BFGS-B", lower = rep(-1,r+m), upper = rep(1,r+m), y = rep(.3,m), D = res$D)
est$par
u1 = rep(0,r+m)
y1 = rep(.3,r+m)
dim(u1)
length(y1)
length(t(D) %*% u1) - y1

## Assume we want to minimize: -(0 5 0) %*% b + 1/2 b^T b
## under the constraints:      A^T b >= b0
## with b0 = (-8,2,0)^T
## and      (-4  2  0) 
##      A = (-3  1 -2)
##          ( 0  0  1)
## we can use solve.QP as follows:
##
Dmat       <- matrix(0,3,3)
diag(Dmat) <- 1
dvec       <- c(0,5,0)
Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
bvec       <- c(-8,2,0)
solve.QP(Dmat,dvec,Amat,bvec=bvec)
