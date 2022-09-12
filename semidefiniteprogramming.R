library(CVXR)

m <- 50

n1 <- global_rho * lambda1
n2 <- global_rho * lambda2
#W <- matrix(.45, m, m)

W <- matrix(runif(m*m), ncol = m)
W <- W %*% t(W) 
W <- W / max(W) 
W <- W 
diag(W) <- 0

D <- create_matrix_D(W, lambda1, lambda2, global_rho)
c <- nrow(D)

-t(D) %*% D

a <- CVXR::Variable(1)

objective <- CVXR::Minimize(a)

constraint1 <- a >= 0
constraint2 <- -t(D) %*% D + diag(a, m) >= 0

problem <- CVXR::Problem(objective, constraints = list(constraint1, constraint2))

solution <- CVXR::solve(problem)

Q = -t(D) %*% D + diag(solution$value + 1, m) 

solution$value
is.positive.semi.definite(Q)
