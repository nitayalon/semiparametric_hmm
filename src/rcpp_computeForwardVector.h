#ifndef rcpp_forwardVector_hpp
#define rcpp_forwardVector_hpp

Rcpp::NumericMatrix rcpp_forwardVector(Rcpp::NumericVector observations, 
                                       Rcpp::NumericMatrix transitionMatrix,
                                 double theta, double p, double b);

#endif