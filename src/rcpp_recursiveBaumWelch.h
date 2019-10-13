#ifndef rcpp_recursiveBaumWelch_hpp
#define rcpp_recursiveBaumWelch_hpp

Rcpp::NumericMatrix rcpp_RecursiveBaumWelch(Rcpp::NumericVector observations,
                                        Rcpp::NumericMatrix transitionMatrix,
                                  double theta, double p, double b);

#endif
