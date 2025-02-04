# Render ----------------------------------------------------------------------------
rmarkdown::render("vignettes/cvn-vignette.Rmd")

# How to build the vignette ---------------------------------------------------------
devtools::build_vignettes()

# How to build the manual ---------------------------------------------------------
devtools::build_manual()




# --- Adjacency matrices
library(CVN)

namen <- paste("m", 1:9, sep="=")
m1 <- create_weight_matrix(type = "grid", 3, 3)
dimnames(m1) <- list(namen, namen)

m2 <- m1
m2[1,5] <- m2[5,1] <- 1
m2[3,5] <- m2[5,3] <- 1
m2[7,5] <- m2[5,7] <- 1
m2[9,5] <- m2[5,9] <- 1
dimnames(m2) <- list(namen, namen)

m3 <- create_weight_matrix(type = "glasso", 3, 3)
dimnames(m3) <- list(namen, namen)

m4 <- create_weight_matrix(type = "full", 2, 3)
dimnames(m4) <- list(namen[1:6], namen[1:6])

save(m1, m2, m3, m4, file = "vignettes/adj.rda")
write.csv2(m1, file = "vignettes/adj1.csv")
write.csv2(m2, file = "vignettes/adj2.csv")
write.csv2(m3, file = "vignettes/adj3.csv")
write.csv2(m4, file = "vignettes/adj4.csv")
# load into excel, save as pdf and make screenshots




