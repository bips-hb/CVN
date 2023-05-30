library(dplyr)
library(CVN)


data(grid)
m <- 9 # must be 9 for this example

#' Choice of the weight matrix W.
#' (uniform random)
W <- matrix(runif(m*m), ncol = m)
W <- W %*% t(W)
W <- W / max(W)
diag(W) <- 0

gamma1 <- 1:10*1e-4
gamma2 <- 1:5*1e-4

lambda1 <- 1:10*1e-4
lambda2 <- 1:5*1e-4

(cvn <- CVN::CVN(grid, W, warmstart = TRUE, eps = 1e-4, maxiter = 1000,
                 gamma1 = gamma1, gamma2 = gamma2, verbose = TRUE))

(cvn$results)

#' if there are less than 3 entries for both lambda1 and lambda2, there is no point
#' in checking whether the optimal estimate based on the AIC or the BIC lies on 
#' the edge of the parameter search space
if(length(gamma1) >= 3 && length(gamma2) >= 3) { 
  
  # write the standard message and fill in the gaps later
  message_warning <- function(criterion = c("BIC", "AIC"), 
                              gamma = c("gamma1", "gamma2"), 
                              hits = c("smallest", "largest")) {
    sprintf("In case you are selecting your model on the basis of the %s: Note that the maximum %s is achieved at the %s value of the %s interval. We suggest you increase the %s interval", 
            criterion[1], 
            criterion[1], 
            hits[1], 
            gamma[1],
            gamma[1])
  }
  
  # For the BIC score ----------------------------------------------------------
  temp <- cvn$results %>% filter(bic == max(bic))
  
  # if gamma1 hits one of the left border
  if (temp$gamma1 == min(gamma1)) {
    warning(message_warning("BIC", "gamma1", "smallest"))
  }
  
  # if gamma1 hits one of the right border
  if (temp$gamma1 == max(gamma1)) {
    warning(message_warning("BIC", "gamma1", "largest"))
  }
  
  # if gamma1 hits one of the left border
  if (temp$gamma2 == min(gamma2)) {
    warning(message_warning("BIC", "gamma2", "largest"))
  }
  
  # if gamma1 hits one of the left border
  if (temp$gamma2 == max(gamma2)) {
    warning(message_warning("BIC", "gamma2", "largest"))
  }
  
  # For the AIC score ----------------------------------------------------------
  temp <- cvn$results %>% filter(aic == max(aic))
  
  # if gamma1 hits one of the left border
  if (temp$gamma1 == min(gamma1)) {
    warning(message_warning("AIC", "gamma1", "smallest"))
  }
  
  # if gamma1 hits one of the right border
  if (temp$gamma1 == max(gamma1)) {
    warning(message_warning("AIC", "gamma1", "largest"))
  }
  
  # if gamma1 hits one of the left border
  if (temp$gamma2 == min(gamma2)) {
    warning(message_warning("AIC", "gamma2", "largest"))
  }
  
  # if gamma1 hits one of the left border
  if (temp$gamma2 == max(gamma2)) {
    warning(message_warning("AIC", "gamma2", "largest"))
  }
}




temp <- cvn$results %>% filter(bic == max(bic))

if (temp$gamma1 == min(gamma1) || temp$gamma2 == max(gamma2)) { 
  
  message(strwrap(
  "In case you are selecting your model on the basis of the BIC:
  Note that the maximum BIC is achieved at the border of the gamma1 
  interval", prefix = " ", initial = ""))
  
  
  warning("WARNING: If you are going by BIC to select your model, the model you are going
          to select ")  
}

if (dplyr::near())

max1 <- gamma1[3]
max2 <- gamma2[3]

expand.grid(
  gamma1 = gamma1, 
  gamma2 = gamma2
)

cat(sprintf("%.1g%s%g\n",min(gamma1),  paste(rep("-", length(gamma1) - 2), collapse = ""), max(gamma1)))

n <- 10
m <- matrix(rbinom(gamma1, 1, .3), ncol = n)


paste(rep(" ", length(gamma1) - 2), collapse = "")
