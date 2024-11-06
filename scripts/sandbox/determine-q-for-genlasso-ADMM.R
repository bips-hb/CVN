library(CVN)
library(matrixcalc)
library(microbenchmark)
library(tictoc)
library(ggplot2)

m = 1000
lambda1 = .5
lambda2 = .4
rho = 1





eta1 = lambda1 / rho
eta2 = lambda2 / rho

W <- CVN::create_weight_matrix(type = "uniform-random", m = m)

DtD = -eta2^2 * (W*W) + eta1^2 * diag(m) + eta2^2*diag(rowSums(W*W)) 

rowsum <- rowSums(W*W)
A = -DtD 

# constructing matrix Q, such that Q - DtD is positive definite
Q = eta1^2*diag(m) + (2*eta2^2)*diag(rowSums(W*W)) 

# get the maximum value on the diag of Q
q = max(eta1^2 + 2*eta2^2 * rowSums(W*W)) 

# determine a with programming
#a <- CVN::matrix_A_inner_ADMM(W, eta1, eta2) 

# there are two versions. 
check <- function(Q, DtD) { 
  A <- Q - DtD
  e <- eigen(A)
  res <- list(
     diag = unique(diag(Q)), 
     Q = Q, 
     A = A, 
     DtD = DtD,
     eigenvalues = e$values, 
     positive_definite = is.positive.definite(A, tol=1e-8)
  )
  class(res) <- "matrixA"
  return(res)
}

print.matrixA <- function(res) { 
  cat(blue("values on the diagonal:\n"))
  for (v in res$diag) { 
    cat(sprintf("%g  ", v))
  }
  cat(blue("\n\neigenvalues: \n"))
  for (e in res$eigenvalues) {
    if (e < 0) {
      cat(red(sprintf("%g  ", e)))
    } else { 
      cat(sprintf("%g  ", e))
    }
  }
  
  if (res$positive_definite) { 
    cat(green(sprintf("\n\n\u2713 positive definite\n")))
  } else { 
    cat(red(sprintf("\n\n\u2717 positive definite\n")))
  }
}

check(Q, DtD) 

check(q*diag(m), DtD)

check(a*diag(m), DtD)

cat(sprintf("\n\tq = %g\ta = %g\n\n", q, a))

e <- eigen(DtD)
C <- max(e$values) * diag(m)
check(C, DtD)

cat(sprintf("\n\tC = %g\ta = %g\n\n", max(e$values), a))

# We can


