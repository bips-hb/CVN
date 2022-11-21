library(CVN)
library(CVNSim)
library(profvis)


p <-  10
m <- 9
lambda1 <- 1
lambda2 <- 3
rho <- 1

eta1 <- lambda1 / rho
eta2 <- lambda2 / rho

set.seed(1)
starting_graph <- CVNSim::generate_graph(p = p, type = "random", probability = .5)
grid_of_graph <- CVNSim::create_grid_of_graphs(starting_graph = starting_graph, 
                                               n_edges_added_x = 2, 
                                               n_edges_removed_x = 1, 
                                               n_edges_added_y = 3, 
                                               n_edges_removed_y = 4)

# generate a single raw dataset for each graph in the grid
set.seed(2)
data <- CVNSim::generate_raw_data_grid(n = 10, grid_of_graph) 

W <- CVNSim::create_weight_matrix("grid")
D <- CVN::create_matrix_D(W, lambda1 = 1, lambda2 = 1, rho = 1)

a <- CVN::matrix_A_inner_ADMM(m, D) + 1

cvn <- CVN(data, W)



Z = updateZ_wrapper(9, p, nrow(D), as.matrix(cvn$Theta[[1]]), as.matrix(cvn$Sigma), W, eta1, eta2, a, 
                rho = 1, maxiter_genlasso = 1000, eps_genlasso = 1e-4, truncate_genlasso = 1e-5)


cvn$Theta[[1]][[1]] + cvn$Sigma[[1]]
Z[[1]] - cvn$Theta[[1]][[1]]

library("Rcpp")
cppFunction('
ListMatrix ListMatrixType(ListMatrix x, ListMatrix y){
            ListMatrix z = x; 
            
            a(0,0) = 100;
            return x;
            }
            ')
x = matrix(list(matrix(0,3,2)),2,2)
a = ListMatrixType(x)
a[[1,1]]
a[[2,2]]