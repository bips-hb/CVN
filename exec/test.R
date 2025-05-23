# Example
fit <- list()
fit$adj_matrices <- vector(length = 2, mode = "list")
fit$adj_matrices[[1]][[1]] <- matrix(c(0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0),
                                nrow = 5, byrow = TRUE)
fit$adj_matrices[[1]][[2]] <- matrix(c(0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0),
                                nrow = 5, byrow = TRUE)
fit$adj_matrices[[1]][[3]] <- matrix(c(0,1,0,1,0,1,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0),
                                nrow = 5, byrow = TRUE)
fit$adj_matrices[[2]][[1]] <- matrix(c(0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0),
                                     nrow = 5, byrow = TRUE)
fit$adj_matrices[[2]][[2]] <- matrix(c(0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0),
                                     nrow = 5, byrow = TRUE)
fit$adj_matrices[[2]][[3]] <- matrix(c(0,1,0,1,1,1,0,1,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0),
                                     nrow = 5, byrow = TRUE)
fit$p <- 5
fit$m <- 3
fit$W <- create_weight_matrix(type="grid", 2,2)
fit$n_lambda_values <- 2
fit$results <- cvn$results[1:2,]
class(fit) <- c("cvn", "list")

out <- visnetwork_cvn(fit, verbose = FALSE)


 core_graphs <- find_core_graph(fit)
 subset_edges <- lapply(core_graphs, function(adj_matrix) {
   as.list(create_edges_visnetwork(adj_matrix))
 })
 
 
 subset_edges <- lapply(core_graphs, function(adj_matrix) {
   create_edges_visnetwork(adj_matrix)
 })
 
 
edges <- create_edges_visnetwork(fit$adj_matrices[[1]][[1]])

edges2 <- set_attributes_to_edges_visnetwork(edges, 
                                            subset_edges = subset_edges[[1]],
                                            width = c(3,1), 
                                            color = c("red", "blue"))


