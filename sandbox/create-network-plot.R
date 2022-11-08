# see https://www.mjdenny.com/Preparing_Network_Data_In_R.html

library(ggraph)
library(igraph)

library(network)
library(sna)
library(ggplot2)


library(GGally)

net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

create_network_object <- function(adj_matrix) { 
   network::as.network(x = as.matrix(adj_matrix), # the network object
                       directed = FALSE, # specify whether the network is directed
                       loops = FALSE, # do we allow self ties (should not allow them)
                       matrix.type = "adjacency") # the type of input
}



M = as.matrix(cvn$adj_matrices[[2]][[1]])

net = create_network_object(cvn$adj_matrices[[2]][[1]])


# vertex names
network.vertex.names(net) = letters[1:10]

set.edge.attribute(net, "core", c(rep(2, 5), rep(.5,5)))
set.edge.attribute(net, "color", c(rep("tomato", 5), rep("black",43-5)))

ggnet2(net, node.size = 12, 
       node.color = "black", 
       edge.size = 1, 
       label = TRUE, 
       label.color = "white", 
       label.size = 5, 
       mode = "circle", 
       edge.color = "color")


bip = data.frame(event1 = c(1, 2, 1, 0),
                 event2 = c(0, 0, 3, 0),
                 event3 = c(1, 1, 0, 4),
                 row.names = letters[1:4])

# weighted bipartite network
bip = network(bip,
              matrix.type = "bipartite",
              ignore.eval = FALSE,
              names.eval = "weights")

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 


ggraph(net) +
  geom_edge_link() +   # add edges to the plot
  geom_node_point()    # add nodes to the plot