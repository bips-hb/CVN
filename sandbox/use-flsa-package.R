library(flsa)
library(CVN)

W <- CVN::create_weight_matrix('grid', m = 9)

W[4, ] <- 0

connList <- vector("list", m)
class(connList) <- "connListObj"

# W <- CVN::create_weight_matrix("full")

lapply(1:m, function(i) {
  indices <- which(W[i, ] != 0)
  if (length(indices) != 0) {
    connList[[i]] <<- as.integer(which(W[i, ] != 0) - 1)
  }
})

names(connList) <- as.character(0:(m-1))

connList_is_null <- all(sapply(connList, function(l) is.null(l)))

connL <- CVN::create_connListObj(W)

if(connL$connList_is_null) {
  fit <- flsa::flsa(1:9, lambda1 = 1, lambda2 = 0)
} else {
  fit <- flsa::flsa(1:9, lambda1 = 1, lambda2 = 1.1, connListObj = connL$connList)
}



# fit <- flsa::flsa(1:9, lambda1 = 1, lambda2 = 1.1, connListObj = connList)
