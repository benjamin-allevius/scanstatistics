#include <Rcpp.h>
using namespace Rcpp;

//' The log probability mass function of the Poisson distribution.
//' 
//' @param x An integer.
//' @param lambda A positive scalar.
//' @return A log-probability.
//' @export
//' @keywords internal
// [[Rcpp::export]]
double poisson_lpmf(const double x, const double lambda) {
  return x * log(lambda) - lgamma(x + 1.0) - lambda;
}
