#include <Rcpp.h>
#include "rcpp_computeCentralDistribution.h"
#include "rcpp_newtonRaphsonTargetFunction.h"
#include "rcpp_computeNewtonRaphsonDenom.h"
using namespace Rcpp;

// [[Rcpp::export]]
double NewtonRaphson(NumericVector data, double theta, double p, int max_iteration)
{
  double b = 0.0;
  
  for(int i = 1; i < max_iteration; ++i)
  {
      double current_target_function_value = optimizationTargetFunction(data,theta,p, b);
      double gradient = rcpp_computeNRDenom(data,theta,p, b);
      double x = current_target_function_value / gradient;
      b = b - x;
  }
  return(b);
}