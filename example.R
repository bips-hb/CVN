####################################
# example.R
#
# Example/structure for running 
# an analysis of a dataset 
####################################

library(CVN)
library(dplyr)
library(readr)

# load the data ---------------------------
data <- CVN::grid 

#' In case that you want to use your own dataset: 
#' data <- readr::read_rds("filename.rds")

# define the weight matrices --------------
W_grid  <- CVN::create_weight_matrix("grid")   # weight matrix for a 3x3 grid
W_empty <- CVN::create_weight_matrix("glasso") # each graph is estimated separately, W = 0 

# function that returns a grid matrix parameterized with alpha
W_radition <- function(alpha) {
  a <- alpha
  b <- 1 - alpha
  matrix(c(
    0, a, 0, 1, 0, 0, 0, 0, 0, 
    a, 0, b, 0, 1, 0, 0, 0, 0, 
    0, b, 0, 0, 0, 1, 0, 0, 0, 
    1, 0, 0, 0, a, 0, 1, 0, 0, 
    0, 1, 0, a, 0, b, 0, 1, 0, 
    0, 0, 1, 0, b, 0, 0, 0, 1, 
    0, 0, 0, 1, 0, 0, 0, a, 0, 
    0, 0, 0, 0, 1, 0, a, 0, b, 
    0, 0, 0, 0, 0, 1, 0, b, 0
  ), 
  ncol = 9)
}

# parameter settings --------------------- 
lambda1 <- 1:3      # sparsity 
lambda2 <- 1:3      # smoothness level

# estimate the CVN -----------------------

# estimate with no smoothing between the graphs
cvn_no_smoothing <- CVN::CVN(data, W = W_empty, 
                             lambda1 = lambda1, lambda2 = lambda2,
                             eps = 1e-4, maxiter = 1e4, verbose = TRUE)

# estimate with no a grid
cvn_grid <- CVN::CVN(data, W = W_grid, 
                     lambda1 = lambda1, lambda2 = lambda2, 
                     eps = 1e-4, maxiter = 1e4, verbose = TRUE)

# estimate with a parameterized graph
cvn_parameterized <- CVN::CVN(data, W = W_radition(0.5), 
                              lambda1 = lambda1, lambda2 = lambda2, 
                              eps = 1e-4, maxiter = 1e4, verbose = TRUE)

# in case you want to plot the results, simply type 
p <- plot(cvn_grid)
p$plots[[1]][[1]]
