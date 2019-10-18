#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <limits>
#include "rcpp_computeAllProbabilities.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_forwardVector(NumericVector observations,
                                 NumericMatrix transitionMatrix,
                                 double theta, double p, double b)
{
  int nObservations = observations.size();
  int nRow = transitionMatrix.nrow();
  int nCol = transitionMatrix.ncol();
  Rcpp::NumericMatrix forward(nObservations,3);

  Rcpp::NumericMatrix pointProbability =
    rcpp_computeAllProbabilities(observations, theta, p, b);
  double initProbability[] = {p * 1.0 /2, 1.0-p, p*1.0/2};

  //init
  for(int i=0; i < 3; ++i)
  {
    forward(0,i) = log(initProbability[i] * pointProbability(0,i));
  }

  for(int k=1;k < nObservations;++k)
  {
    for(int state=0;state<nRow; ++state)
    {
      double logsum = -std::numeric_limits<double>::infinity();
      for(int previousState=0; previousState < nCol; ++previousState)
      {
        double temp   = forward(k-1,previousState) +
          log(transitionMatrix(previousState,state));
        if(temp > -std::numeric_limits<double>::infinity())
        {
          logsum = temp + log(1 + exp(logsum - temp ));
        }
      }
      forward(k,state) = log(pointProbability(k,state)) + logsum;
    }
  }
  return forward;
}
