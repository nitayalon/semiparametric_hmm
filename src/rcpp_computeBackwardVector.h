#ifndef rcpp_backwardVector_hpp
#define rcpp_backwardVector_hpp

Rcpp::NumericMatrix rcpp_backwardVector(Rcpp::NumericVector observations, 
                                        Rcpp::NumericMatrix transitionMatrix,
                                  double theta, double p, double b);
  
#endif