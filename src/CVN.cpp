#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::DoubleVector aug_genlassoRcpp(Rcpp::NumericVector y, 
                                const Rcpp::NumericMatrix W, 
                                const int m, 
                                const int c, 
                                const double lambda1, 
                                const double lambda2, 
                                const double global_rho, 
                                double a, 
                                const double rho, 
                                const int max_iter,
                                const double eps) { 
  int i,j,k ;
  
  // initialize vectors
  Rcpp::DoubleVector beta_new (m) ; 
  Rcpp::DoubleVector beta_old (m) ; 
  
  Rcpp::DoubleVector alpha_new (c) ; 
  Rcpp::DoubleVector alpha_old1 (c) ;
  Rcpp::DoubleVector alpha_old2 (c) ;
  Rcpp::DoubleVector alpha (c) ;
  
  Rcpp::DoubleVector delta (m) ;
  
  // initialize some values that are used repeatedly 
  double eta1 = global_rho * lambda1 ; 
  double eta2 = global_rho * lambda2 ; 
  a = rho*a ; 
  double C = 1 / (1 + a) ; 
  
  //Rprintf("%f, %f, %f, %f\n", eta1, eta2, a, C) ; 
  
  Rcpp::DoubleVector ya = C * y ; 
  
  int iter = 0 ; 
  double diff = 0 ;
  
  while (iter < max_iter) { 
    Rprintf("iter: %d\n", iter) ; 
    
    alpha = 2*alpha_old1 - alpha_old2 ; 
    
    for(i = 0; i < m; i++) { 
      delta[i] = eta1 * alpha[i] ;  
      
      for (j = i; j < m; j++) { 
         delta[i] = delta[i] + eta2*W[i,j]*alpha[i + j + m - 1] ; 
      }
      
      for (j = 0; j < i; j++) { 
        delta[i] = delta[i] - eta2*W[j,i]*alpha[i + j + m - 1] ; 
      }
      Rprintf("delta[%d]: %f\n", i, delta[i]) ; 
    }
    
    beta_new = (1 - C)*beta_old + ya - C*delta ; 
    
    for (i = 0; i < m; i ++) { 
      Rprintf("%f\t%f\t%f\n", beta_new[i], beta_old[i], ya[i]) ;  
    }
    
    for (i = 0; i < m; i ++) { 
      diff += abs(beta_new[i] - beta_old[i]) ;  
    }
    
    if (diff < eps) { 
      Rprintf("diff = %f\teps: %f\n", diff, eps) ; 
      return(beta_new) ;  
    }
    
    diff = 0 ; 
    
    for (i = 0; i < c; i++) { 
      alpha_new[i] = alpha_old1[i] ;   
    }
    
    for (i = 0; i < m; i++) { 
       alpha_new[i] = eta1*beta_new[i] ; 
    }
    
    k = m; 
    for (i = 0; i < m-1; i++) { 
      for (j = i+1; j < m; j++) { 
        alpha_new[k] = rho * eta2*W[i,j]*(beta_new[i] - beta_new[j]) ; 
        if (alpha_new[k] > 1) { 
          alpha_new[k] = 1 ;  
        } 
        if (alpha_new[k] < -1) { 
          alpha_new[k] = -1 ;  
        } 
        k++; 
      }
    }
//alpha_new <- alpha_old1 + rho * D %*% beta_new
    
    for (i = 0; i < m; i++) { 
      beta_old[i]   = beta_new[i] ; 
      alpha_old2[i] = alpha_old1[i] ; 
      alpha_old1[i] = alpha_new[i] ; 
    }
    
    
    iter ++; 
  }
  
  return(beta_new) ;   
}

// [[Rcpp::export]]
double fisherTestGreater(int a, int b, int c, int d) {

  // compute the marginal counts for the drug and the event
  int drug_marg  = a + c ;
  int event_marg = a + b ;
  int n          = a + b + c + d ;

  // probability of observing the actual table
  double p_obs_table = R::dhyper(a, drug_marg, n - drug_marg, event_marg, false) ;

  double p_value = 0.0 ;
  double p_table ;

  // walk through all possible tables
  int max_a = drug_marg ;
  if (max_a > event_marg) {
    max_a = event_marg ;
  }

  for (int i = a; i <= max_a; i ++) {
    p_table = R::dhyper(i, drug_marg, n - drug_marg, event_marg, false) ;
    if (p_table <= p_obs_table) {
      p_value += p_table ;
    }
  }

  return p_value ;
}

// [[Rcpp::export]]
double midPFisherTestGreater(int a, int b, int c, int d) {

  // compute the marginal counts for the drug and the event
  int drug_marg  = a + c ;
  int event_marg = a + b ;
  int n          = a + b + c + d ;

  // probability of observing the actual table
  double p_obs_table = R::dhyper(a, drug_marg, n - drug_marg, event_marg, false) ;
  return fisherTestGreater(a, b, c, d) - p_obs_table / 2.0 ;
}


//' Create 2 x 2 Tables
//'
//' Creates a data frame containing all 2 x 2 contingency tables
//' given a raw spontaneous reporting (SR) data set. An SR data set
//' is a binary matrix, where each row is a report. The first
//' columns represent the presence or absence of a drug, the
//' The other columns represent the presence or absence of an event.
//' See for more information the wrapper function,
//' \code{\link{convertRawReports2Tables}}.
//'
//' The code is a simplified version of the function \code{create2x2TablesRcpp}
//' in the \code{SRSim} package.
//'
//' @param reports A binary matrix. Each row is a report
//' @param n_drugs The number of drugs
//' @param n_events The number of events
//'
//' @return A dataframe. A description of the columns can be found in the commentary
//'         for the function \code{\link{convertRawReports2Tables}}
//'
//' @seealso \code{\link{convertRawReports2Tables}}
// [[Rcpp::export]]
Rcpp::DataFrame convertRawReports2TablesRcpp (Rcpp::IntegerMatrix reports, int n_drugs, int n_events) {

  int n_pairs = n_drugs * n_events ;
  int k ;

  // create vectors that will make up the data frame
  Rcpp::IntegerVector drug_id (n_pairs) ;
  Rcpp::IntegerVector event_id (n_pairs) ;
  Rcpp::IntegerVector a (n_pairs, 0) ;
  Rcpp::IntegerVector b (n_pairs, 0) ;
  Rcpp::IntegerVector c (n_pairs, 0) ;
  Rcpp::IntegerVector d (n_pairs, 0) ;

  // run over all pairs and fill the vectors
  for (int i = 0; i < n_drugs; i ++) {
    for (int j = 0; j < n_events; j ++) {
      k = i*n_events + j ; // current pair index
      drug_id[k]        = i+1 ;
      event_id[k]       = j+1 ;
    }
  }

  // compute the tables
  // go over all the reports
  bool drug, event ;
  for (int r = 0; r < reports.nrow(); r ++) {
    // go over all drug-event pairs
    for (int i = 0; i < n_drugs; i ++) {
      for (int j = 0; j < n_events; j ++) {
        k = i*n_events + j ; // pair index
        drug = (reports(r, i) == 1) ;
        event = (reports(r, n_drugs + j) == 1) ;
        if (drug) {
          if (event) {
            a[k] ++ ;
          } else {
            c[k] ++ ;
          }
        } else {
          if (event) {
            b[k] ++ ;
          } else {
            d[k] ++ ;
          }
        }
      }
    }
  }

  return Rcpp::DataFrame::create( Named("drug_id") = drug_id,
                                  Named("event_id") = event_id,
                                  Named("a") = a,
                                  Named("b") = b,
                                  Named("c") = c,
                                  Named("d") = d
  ) ;
}
