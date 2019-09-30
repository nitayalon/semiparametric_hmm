#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <limits>
#include "rcpp_computeAllProbabilities.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_backwardVector(NumericVector observations, 
                                 NumericMatrix transitionMatrix,
                                 double theta, double p, double b)
{
  int nObservations = observations.size();
  int nRow = transitionMatrix.nrow();
  int nCol = transitionMatrix.ncol();
  Rcpp::NumericMatrix backword(nObservations,3);
  
  Rcpp::NumericMatrix pointProbability = 
    rcpp_computeAllProbabilities(observations, theta, p, b);
  //init
  for(int i=0; i < 3; ++i)
  {
    backword(nObservations-1,i) = 0.0;
  }
  int m = nObservations-1;
  
  for(int k=m-1;k >= 0;--k)
  {
    for(int state = 0;state < nRow; ++state)
    {
      double logsum = -std::numeric_limits<double>::infinity();
      
      for(int nextState = 0; nextState < nCol; ++nextState)
      {
        double x = backword(k+1,nextState);
        
        double temp  = backword(k+1,nextState) +
          log(transitionMatrix(state, nextState) * pointProbability(k+1,state));

        if(temp > -std::numeric_limits<double>::infinity())
        {
          logsum = temp + log(1 + exp(logsum - temp));
        }
      }
      backword(k,state) = logsum;
    }
  }
  return backword;
}
