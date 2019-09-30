#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double computeCentralDistribution(double x, double theta, double p, double b) {
  double g_x = 1 / (1 - p + p * (exp(-1*b) * cosh(theta * x)));
  return g_x;
}
