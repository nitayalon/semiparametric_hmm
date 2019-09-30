#include <Rcpp.h>
#include "rcpp_computePointDensity.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_computeAllProbabilities(NumericVector observations, 
                                           double theta, 
                                           double p, 
                                           double b)
{
  int N = observations.size();
  Rcpp::NumericMatrix probabilityMatrix(N,3);
  Rcpp::NumericVector driftParameter(3);
  Rcpp::NumericVector pointProbability(N);
  
  for(int k = 0; k < 3; ++k)
  {
    driftParameter[k] = (k-1) * theta;
  }
  
  for(int i = 0; i < N; ++i)
  {
    pointProbability(i) = rcpp_pointDensity(observations(i),
                     theta, p, b, N);
  }
  
  for(int i = 0; i < N; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      probabilityMatrix(i,j) = pointProbability(i) * 
        exp(driftParameter(j) * observations(i));
    }
  }
  return probabilityMatrix;
}
