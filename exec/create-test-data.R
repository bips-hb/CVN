#' Requires the installation of CVNSim, see https://github.com/bips-hb/CVNSim

library(CVNSim)

p <- 10  # number of variables
m <- 9   # number of graphs. !!! must be 9 !!!

# Choice of the weight matrix W
uniform_random <- F   
grid_weight <- T
predefined <- F 
disconnected <- F

W <- matrix(1, m, m) # standard fully-connected 

if (disconnected) { 
  W <- matrix(0, m, m)
}

if (uniform_random) { 
  W <- matrix(runif(m*m), ncol = m)
  W <- W %*% t(W)
  W <- W / max(W)
  W <- W
  diag(W) <- 0
}

if (grid_weight) { 
  W <- CVNSim::create_weight_matrix(type = "grid") 
}

if (predefined) { 
  e <- 0
  W <- matrix(c(0, 1, e, 1,
                1, 0, 1, e,
                e, 1, 0, 1,
                1, e, 1, 0), ncol = 4)
}


starting_graph <- CVNSim::generate_graph(p, type = "random", probability = .5)
grid3 <- create_grid_of_graphs(starting_graph = starting_graph, 
                               n_edges_added_x = 2, 
                               n_edges_removed_x = 2, 
                               n_edges_added_y = 2,
                               n_edges_removed_y = 2,
                               verbose = TRUE)

grid <- CVNSim::generate_raw_data_grid(100, grid3)
save(grid, file="data/grid.RData")
