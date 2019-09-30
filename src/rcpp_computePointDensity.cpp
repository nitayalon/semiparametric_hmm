#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double rcpp_pointDensity(double x, double theta, double p, double b, int N) {
  double point_density = (1 / (N * 1.0)) * 
    (1 / (1 - p + p * exp(-1 * b * cosh(theta * x))));
  return point_density;
}

