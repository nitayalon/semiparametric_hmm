#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <limits>
#include "rcpp_computeAllProbabilities.h"
#include "rcpp_computeForwardVector.h"
#include "rcpp_computeBackwardVector.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_RecursiveBaumWelch(NumericVector observations,
                             NumericMatrix transitionMatrix,
                             double theta, double p, double b)
{
  int nObservations = observations.size();
  int nRow = transitionMatrix.nrow();
  int nCol = transitionMatrix.ncol();

  Rcpp::NumericMatrix pointProbability =
    rcpp_computeAllProbabilities(observations, theta, p, b);
  Rcpp::NumericMatrix forward_vector = rcpp_forwardVector(observations, transitionMatrix,
                                                          theta, p, b);
  Rcpp::NumericMatrix backward_vector = rcpp_backwardVector(observations, transitionMatrix,
                                                    theta, p, b);
  double probObservations = forward_vector((nObservations-1),0);

  for(int i=1;i<nRow;++i)
  {
    double j = forward_vector((nObservations-1),i);
    if(j > -std::numeric_limits<double>::infinity())
    {
      probObservations = j + log(1+exp(probObservations-j));
    }
  }

  for(int x=0; x<nRow; ++x)
  {
    for(int y=0; y<nRow; ++y)
    {
      double temp = forward_vector(x,0) + log(transitionMatrix(x,y)) +
        log(pointProbability(y,1)) + backward_vector(y,1);

      for(int i=1;i<(nObservations-2);++i)
      {
        double j = forward_vector(i,x) + log(transitionMatrix(x,y)) +
          log(pointProbability(i+1,y)) + backward_vector(i+1,y);
        if(j > -std::numeric_limits<double>::infinity())
        {
          temp = j + log(1+exp(temp-j));
        }
      }

      temp = exp(temp - probObservations);
      transitionMatrix(x,y) = temp;
    }
  }
  return(transitionMatrix);
}


