# Render ----------------------------------------------------------------------------
rmarkdown::render("vignettes/cvn-vignette.Rmd")

# How to build the vignette ---------------------------------------------------------
devtools::build_vignettes()

# How to build the manual ---------------------------------------------------------
devtools::build_manual()




# cvnfit.RData ----------------------------------------------------------------------
library(CVN)
library(dplyr)

# Simulate the dataset
set.seed(2024)  
n <- 300  

# Create 10 normally distributed variables for the graph
data <- as.data.frame(matrix(rnorm(n * 10), ncol = 10))
colnames(data) <- paste0("X", 1:10)

# Add two discrete external covariates
data$income <- sample(c("low","middle", "high"), n, replace = TRUE)  
data$bmi <- sample(c("underweight","normal weight", "overweight"), n, replace = TRUE)

# Split the dataset into subsets based on dosis and bmi
data_list <- data %>%
  group_by(income, bmi) %>%
  group_split() %>%
  lapply(function(df) df %>% select(-income, -bmi))

names(data_list) <- 
  apply(expand.grid(income = c("low","middle", "high"), 
                    bmi = c("underweight","normal weight", "overweight")), 
        1,
        function(x) paste0(x[1], "_", x[2]))

W_grid <- create_weight_matrix(type = "grid", k = 3, l = 3, plot = FALSE)

lambda1 = seq(0.5, 2, length = 3)  # sparsity
lambda2 = c(1, 1.5)                # smoothing

cvn <- CVN(data = data_list, 
           W = W_grid, 
           lambda1 = lambda1, 
           lambda2 = lambda2, 
           eps = 1e-2,       # makes it faster but less precise; default = 1e-4
           maxiter = 500, 
           n_cores = 1,          # no parallizing
           warmstart = TRUE,       # uses the glasso
           verbose = FALSE)

fit <- extract_cvn(cvn, 6)
save(fit, file = "inst/cvnfit.RData")

# interpolate
interpolate6 <- interpolate(fit6, c(0,0,0,0,0,0,0,0.5,0.5), truncate = 0.05)
fit10 <- combine_cvn_interpolated(fit6, interpolate6)

fit10 <- visnetwork_cvn(fit10)













# Adjacency matrices for the vignette ------------------------------------------
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




