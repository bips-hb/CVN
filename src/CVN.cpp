#include <Rcpp.h>
using namespace Rcpp;

//' Solving Generalized LASSO with fixed \eqn{\lambda = 1}
//' 
//' Solves efficiently the generalized LASSO problem of the form 
//' \deqn{
//'   \hat{\beta} = \text{argmin } \frac{1}{2} || y - \beta ||_2^2 + ||D\beta||_1 
//' }
//' where \eqn{\beta} and \eqn{y} are \eqn{m}-dimensional vectors and 
//' \eqn{D} is a \eqn{(c \times m)}-matrix where \eqn{c \geq m}. 
//' We solve this optimization problem using an adaption of the ADMM
//' algorithm presented in Zhu (2017). 
//'
//' @return The estimated vector \eqn{\hat{\beta}}
//'
//' @references 
//' Zhu, Y. (2017). An Augmented ADMM Algorithm With Application to the 
//' Generalized Lasso Problem. Journal of Computational and Graphical Statistics, 
//' 26(1), 195–204. https://doi.org/10.1080/10618600.2015.1114491
//' 
//' @seealso \code{\link{convertRawReports2Tables}}
// [[Rcpp::export]]
Rcpp::NumericVector aug_genlassoRcpp(Rcpp::NumericVector y, 
                                const Rcpp::NumericMatrix& W, 
                                const int m, 
                                const int c, 
                                const double eta1, 
                                const double eta2, 
                                double a, 
                                const double rho, 
                                const int max_iter,
                                const double eps) { 
  int i,j,k ; // indices
  
  /* some frequently used constants */
  a = rho*a ; 
  double C = 1 / (1 + a) ; 
  
  //Rf_PrintValue(W) ; 
  
  /* initialize vectors for beta-update step in the ADMM */
  Rcpp::NumericVector beta_new (m) ; // beta^(k + 1)
  Rcpp::NumericVector beta_old (m) ; // beta^k
  Rcpp::NumericVector delta (m) ; // aux. vector for beta-update step

  /* initialize vectors for alpha-update step in the ADMM */
  Rcpp::NumericVector alpha_new (c) ; 
  Rcpp::NumericVector alpha_old1 (c) ;
  Rcpp::NumericVector alpha_old2 (c) ;
  Rcpp::NumericVector alpha (c) ; // used to store (2*alpha^k - alpha^(k-1))

  /* Indices used for the alpha-update step */
  Rcpp::IntegerVector steps (m-1) ; 
  steps[0] = 0 ; 
  
  for (i = 1; i < (m-1); i ++) { 
    steps[i] = m - i + steps[i - 1] ; 
  }
  
  int iter = 0 ;    // number of iterations 
  double diff = 0 ; // absolute difference between beta^(k+1) and beta^k
  
  /* loop until either the max. no. of iterations are reached 
   * or the difference (diff) is smaller then eps
   */
  while (iter < max_iter) { 
    
    /* ------- beta-update step ---------*/
    for (i = 0; i < c; i++) { 
      alpha[i] = 2*alpha_old1[i] - alpha_old2[i] ; 
    }
    
    // go over all possible pairs (i,j), same as D^T %*% (2 alpha^(k) - alpha^(k-1)) 
    for (i = 0; i < m; i++) { 
      delta[i] = eta1 * alpha[i] ;  
      
      for (j = i+1; j < m; j++) { 
         delta[i] = delta[i] + eta2*W(i,j)*alpha[m + steps[i] + (j - i) - 1] ; 
      }
      
      for (j = 0; j < i; j++) { 
        delta[i] = delta[i] - eta2*W(j,i)*alpha[m + steps[j] - (j - i) - 1] ; 
      }
    }
    
    //Rf_PrintValue(delta);
    
    // update beta with the computed delta and determine difference
    diff = 0 ;
    for (i = 0; i < m; i ++) { 
      beta_new[i] = C*(a*beta_old[i] + y[i] - delta[i]) ;
      diff += abs(beta_new[i] - beta_old[i]) ;  
    }
    
    // determine whether converged or not
    if (diff < eps) { 
      /* Turn to zero when really close */
      for (i = 0; i < m; i ++) { 
        if (abs(beta_new[i]) < 10e-7) { 
          beta_new[i] = 0 ;  
        } 
      }
      return(beta_new) ;  
    }
    
    /* --------- alpha update step ----------- */
    for (i = 0; i < m; i++) { 
       alpha_new[i] = alpha_old1[i] + rho * eta1 * beta_new[i] ; 
    }
    
    k = m; 
    // go over all unique pairs (i,j)
    for (i = 0; i < m-1; i++) { 
      for (j = i+1; j < m; j++) { 
        alpha_new[k] = alpha_old1[k] + rho * eta2*W(i,j)*(beta_new[i] - beta_new[j]) ; 
        //Rprintf("k = %d\t(i,j) = (%d,%d)\tW[%d,%d] = %g\t%g --> %g\n", k, i, j, i+1, j+1, W[i,j], alpha_old1[k], alpha_new[k]) ; 
        k++; 
      }
    }
    
    /* Threshold alpha. Must lie in [-1, 1] range */ 
    for (i = 0; i < c; i ++) { 
      if (alpha_new[i] > 1) { 
        alpha_new[i] = 1 ;  
      } 
      if (alpha_new[i] < -1) { 
        alpha_new[i] = -1 ;  
      } 
    }
    
    /* update beta and alpha for the next iteration step */
    beta_old = Rcpp::clone(beta_new) ; 
    alpha_old2 = Rcpp::clone(alpha_old1) ; 
    alpha_old1 = Rcpp::clone(alpha_new) ; 
    
    iter ++; 
  }
  
  /* Turn to zero when really close */
  for (i = 0; i < m; i ++) { 
    if (abs(beta_new[i]) < 10e-7) { 
      beta_new[i] = 0 ;  
    } 
  }
  
  return(beta_new) ;   
}
