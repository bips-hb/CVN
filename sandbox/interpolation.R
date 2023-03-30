source("exec/quick-testrun.R")

#weights <- runif(9)
weights <- rep(1, 9)
lambda1 = 1:2
lambda2 = 1:2 
gamma1 = NULL
gamma2 = NULL
truncate = NULL

cvn$Theta

lambda1 <- cvn$results$lambda1
lambda2 <- cvn$results$lambda2

# Function for interpolation
inter <- function(x, y, w, lambda1, lambda2) { 
    lambda1 * abs(x) + lambda2 * sum(abs(w*y - x)) 
}

# Go over all lambda1 and lambda2 pairs
lapply(1:cvn$n_lambda_values, function(i) {
  lambda1 <- cvn$results$lambda1[i]
  lambda2 <- cvn$results$lambda1[i]
  
  # initial the matrix
  adj_matrix <- Matrix(0, nrow = cvn$p, ncol = cvn$p, sparse = TRUE)
  
  # go over all potential edges (s,t)
  combn(cvn$p, 2, function(pair) {
    y <- sapply(cvn$Theta[[1]], function(Theta_i) Theta_i[pair[1],pair[2]])
    
    fit <- optim(
      par = 0,
      fn = inter,
      y = y,
      w = weights,
      lambda1 = lambda1,
      lambda2 = lambda2,
      method = "Brent",
      lower = min(weights * y) - 1,
      upper = max(weights * y) + 1
    )
    
    edge_exists <- abs(fit$par) >= truncate
    cat(sprintf("par: %g\tvalue:%g\t--->\t%d\n", fit$par, fit$value, abs(fit$par) >= truncate))
    adj_matrix[pair[1], pair[2]] <<- as.numeric(edge_exists)
    adj_matrix[pair[2], pair[1]] <<- as.numeric(edge_exists)
  })
  
  return(adj_matrix)
})

combn(4, 2, function(i) i + j)

int = CVN::interpolate(cvn, weights = rep(1,9))




library(CVNSim)

set.seed(1)
starting_graph <- CVNSim::generate_graph(p = 20, type = "random", probability = .1)
grid_of_graph <- CVNSim::create_grid_of_graphs(starting_graph = starting_graph, 
                                               n_edges_added_x = 2, 
                                               n_edges_removed_x = 1, 
                                               n_edges_added_y = 3, 
                                               n_edges_removed_y = 4)

# generate a single raw dataset for each graph in the grid
set.seed(2)
data <- CVNSim::generate_raw_data_grid(n = 100, grid_of_graph) 

cvn <- CVN::CVN(data, W = CVN::create_weight_matrix("full"))

int <- interpolate(cvn, weights = rep(100,9))
cvn$adj_matrices[[3]][[9]]
int[[3]]

differ <- function(adj) { 
  combn(4,2, function(pair) {
    sum(abs(adj[[pair[1]]] - adj[[pair[2]]]))
  })
}

differ(int)
differ(cvn$adj_matrices[[4]])

combn(4,2, function(pair) {
  sum(abs(int[[pair[1]]] - int[[pair[2]]]))
})

combn(4,2, function(pair) {
  sum(abs(int[[pair[1]]] - int[[pair[2]]]))
})

y <- sapply(cvn$Theta[[1]], function(Theta_i) Theta_i[1,2])

if (sum(y) == 0) { # no edge
  cat("no edge\n")
}


fit <- optim(
  par = 0,
  fn = inter,
  y = y,
  w = weights,
  lambda1 = lambda1[1],
  lambda2 = lambda2[1],
  method = "Brent",
  lower = min(weights * y) - 1,
  upper = max(weights * y) + 1
)




