#include <numeric>
#include <Rcpp.h>
#include "rcpp_computeCentralDistribution.h"
using namespace Rcpp;

// [[Rcpp::export]]
double optimizationTargetFunction(NumericVector x, double theta, double p, double b) {
  int N = x.size();
  NumericVector target_function(N);
  for(int i = 0; i < N; ++i) {
    target_function[i] = computeCentralDistribution(x[i], theta, p, b);
  }
  double target_function_value = std::accumulate(target_function.begin(), 
                                                 target_function.end(), 0.0);
  return(target_function_value - N);
}