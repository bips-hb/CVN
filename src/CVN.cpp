#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector aug_genlassoRcpp(Rcpp::NumericVector y, 
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
  Rcpp::NumericVector beta_new (m) ; 
  Rcpp::NumericVector beta_old (m) ; 
  
  Rcpp::NumericVector ya (m) ; 
  
  Rcpp::NumericVector alpha_new (c) ; 
  Rcpp::NumericVector alpha_old1 (c) ;
  Rcpp::NumericVector alpha_old2 (c) ;
  Rcpp::NumericVector alpha (c) ;
  
  for (i = 1; i < m; i ++) { 
    alpha_old1[i] = 0 ;  
  }
  
  Rcpp::NumericVector delta (m) ;
  
  // initialize some values that are used repeatedly 
  Rcpp::IntegerVector steps (m) ; 
  steps[0] = 0 ; 
  
  for (i = 1; i < m; i ++) { 
    steps[i] = m - i + steps[i - 1] ; 
    //Rprintf("%d ", steps[i]) ; 
  }
  
  double eta1 = lambda1 / global_rho ; 
  double eta2 = lambda2 / global_rho ; 
  a = rho*a ; 
  double C = 1 / (1 + a) ; 
  
  Rprintf("%f, %f, %f, %f\n", eta1, eta2, a, C) ; 
  
  for (i = 0; i < m; i ++) { 
     ya[i] = C*y[i] ; 
  }
  
  int iter = 0 ; 
  double diff = 0 ;
  
  while (iter < max_iter) { 
    //Rprintf("iter: %d\n", iter) ; 
    
    for (i = 0; i < c; i++) { 
      alpha[i] = 2*alpha_old1[i] - alpha_old2[i] ; 
    }
    
    for(i = 0; i < m; i++) { 
      //Rprintf("i = %d -----------\n", i) ; 
      delta[i] = eta1 * alpha[i] ;  
      
      for (j = i+1; j < m; j++) { 
         //Rprintf("+ (%d, %d)\t%d\n", i, j, m + steps[i] + (j - i) - 1) ;   
         delta[i] = delta[i] + eta2*W[i,j]*alpha[m + steps[i] + (j - i) - 1] ; 
      }
      
      for (j = 0; j < i; j++) { 
        //Rprintf("- (%d, %d)\t%d\n", i, j, m + steps[j] - (j - i) - 1) ;   
        delta[i] = delta[i] - eta2*W[j,i]*alpha[m + steps[j] - (j - i) - 1] ; 
      }
      
      delta[i] = C*delta[i] ; 
      //Rprintf("delta[%d]: %f\n", i, delta[i]) ; 
    }
    
    for (i = 0; i < m; i ++) { 
      beta_new[i] = a*beta_old[i] + ya[i] - delta[i] ; 
      //Rprintf("%g  ", beta_new[i]) ; 
    }
    //Rprintf("\n") ; 
    //for (i = 0; i < m; i ++) { 
    //  Rprintf("%f\t%f\t%f\n", beta_new[i], beta_old[i], ya[i]) ;  
    //}
    //Rf_PrintValue(beta_new);
    //Rf_PrintValue(beta_old);
    
    diff = 0 ; 
    for (i = 0; i < m; i ++) { 
      diff += abs(beta_new[i] - beta_old[i]) ;  
    }
    
    Rprintf("diff: %g\n", diff) ; 
    
    //if(iter == 2) { 
    //  return(beta_new) ;  
    //}
    
    if (diff < eps) { 
      Rprintf("finished. iter = %d\tdiff = %f\teps: %f\n", iter, diff, eps) ; 
      return(beta_new) ;  
    }
    
    /*for (i = 0; i < c; i++) { 
      alpha_new[i] = alpha_old1[i] ;   
    }*/
    
    for (i = 0; i < m; i++) { 
       alpha_new[i] = alpha_old1[i] + rho * eta1 * beta_new[i] ; 
    }
    
    k = m; 
    for (i = 0; i < m-1; i++) { 
      for (j = i+1; j < m; j++) { 
        //Rprintf("k: %d\t (i,j) = (%d,%d)\n", k, i, j) ;
        alpha_new[k] = alpha_old1[k] + rho * eta2*W[i,j]*(beta_new[i] - beta_new[j]) ; 
        if (alpha_new[k] > 1) { 
          alpha_new[k] = 1 ;  
        } 
        if (alpha_new[k] < -1) { 
          alpha_new[k] = -1 ;  
        } 
        k++; 
      }
    }
    
    Rf_PrintValue(alpha_new) ; 
    Rf_PrintValue(alpha_old1) ;
    Rf_PrintValue(alpha_old2) ;
    Rprintf("\n") ; 
    
//alpha_new <- alpha_old1 + rho * D %*% beta_new
    
    /*for (i = 0; i < m; i++) { 
      beta_old[i]   = beta_new[i] ; 
    }*/
    
    //std::copy( beta_new.begin(), beta_new.end(), beta_old.begin() ) ;
    //std::copy( alpha_old2.begin(), alpha_old2.end(), alpha_old1.begin() ) ;
    //std::copy( alpha_new.begin(), alpha_new.end(), alpha_old2.begin() ) ;
    
    /*Rprintf("\n\nbefore:\n") ; 
    Rf_PrintValue(beta_new);
    Rf_PrintValue(beta_old);
    Rf_PrintValue(alpha_new);
    Rf_PrintValue(alpha_old1);
    Rf_PrintValue(alpha_old2);*/
    
    /*for (i = 0; i < m; i ++) { 
      beta_old[i] = beta_new[i] ;  
    }
    
    for (i = 0; i < c; i ++) { 
      alpha_old2[i] = alpha_old1[i] ; 
      alpha_old1[i] = alpha_new[i] ; 
    }*/
    
    beta_old = Rcpp::clone(beta_new) ; 
    alpha_old2 = Rcpp::clone(alpha_old1) ; 
    alpha_old1 = Rcpp::clone(alpha_new) ; 
    
    /*Rprintf("\n\nafter:\n") ; 
    Rf_PrintValue(beta_new);
    Rf_PrintValue(beta_old);
    Rf_PrintValue(alpha_new);
    Rf_PrintValue(alpha_old1);
    Rf_PrintValue(alpha_old2);*?
    
    /*for (i = 0; i < c; i++) { 
      alpha_old2[i] = alpha_old1[i] ; 
      alpha_old1[i] = alpha_new[i] ; 
    }
    
    Rprintf("beta_new:\n") ; 
    for(i = 0; i < m; i ++) { 
      Rprintf("%g  ", beta_new[i]) ; 
    }
    Rprintf("\n") ; */
    
    if (iter == 2) { 
      return(beta_old) ; 
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
