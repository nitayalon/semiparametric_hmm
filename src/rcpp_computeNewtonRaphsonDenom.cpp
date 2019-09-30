#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double rcpp_computeNRDenom(NumericVector x, double theta, double p, double b)
{
  int N = x.size();
  NumericVector numerator = p * exp(-b) * cosh(theta*x);
  NumericVector denom = pow((1 - p + p * exp(-b) * cosh(theta*x)),2.0);
  NumericVector nr_denom(N);
  for (int i = 0; i < N; ++i) {
    nr_denom[i] =   numerator[i] / denom[i];
  }
  double nr_denom_sum = std::accumulate(nr_denom.begin(), nr_denom.end(), 0.0);
  return(nr_denom_sum);
}