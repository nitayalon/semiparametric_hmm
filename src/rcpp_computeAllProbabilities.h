#ifndef rcpp_computeAllProbabilities_hpp
#define rcpp_computeAllProbabilities_hpp

Rcpp::NumericMatrix rcpp_computeAllProbabilities(Rcpp::NumericVector observations, 
                                           double theta, 
                                           double p, 
                                           double b);

#endif