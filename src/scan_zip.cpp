#include <Rcpp.h>
#include <cmath>
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

//' Terms in the log-product (sum) of the EB-ZIP window statistic.
//' 
//' Computes one or more terms in the log-product (sum) of the numerator or 
//' denominator of the EB-ZIP window statistic (i.e. the one calculated for a 
//' given space-time window W).
//' @param p Numeric vector of structural zero probabilities.
//' @param d Numeric vector of estimates of the structural zero indicators. Of 
//'    same length as \code{p}.
//' @param mu Numeric vector of given/estimated Poisson expected value 
//' parameters. Of same length as \code{p}.
//' @param y Integer vector of observed counts. Of same length as \code{p}.
//' @param tol Scalar; probability p below this is considered equal to zero.
//' @return A numeric vector of same length as input vector \code{p}.
//' @keywords internal
// [[Rcpp::export]]
NumericVector zip_statistic_logfactor(const NumericVector& p, 
                                      const NumericVector& d,
                                      const NumericVector& mu,
                                      const NumericVector& y,
                                      double tol = 1e-08) {
  int n = y.size();
  NumericVector res(n);
  for (int i = 0; i < n; ++i) {
    if (p[i] < tol) {
      res[i] = poisson_lpmf(y[i], mu[i]);
    } else {
      res[i] = d[i] * log(p[i]) + (1 - d[i]) * (log(1 - p[i]) + poisson_lpmf(y[i], mu[i]));
    }
  }
  return res;
}

//' Estimate the relative risk for an outbreak using a ZIP distribution.
//' 
//' Scalar estimate of the relative risk \eqn{q} for zero-inflated Poisson data.
//' @param d A vector of indicator variables for whether the corresponding count
//'    in the argument \code{y} is an excess zero or not. Can also be estimates 
//'    of these indicators.
//' @param mu A vector of given/estimated Poisson expected value parameters, of 
//'    same length as \code{d}.
//' @param y An integer vector of the observed counts, of same length as 
//'    \code{d}.
//' @return A scalar, the estimated relative risk.
//' @keywords internal
// [[Rcpp::export]]
double estimate_zip_relrisk(const NumericVector& d,
                            const NumericVector& mu,
                            const NumericVector& y) {
  int n = y.size();
  double q;
  double numerator = 0.0;
  double denominator = 0.0;
  
  for (int i = 0; i < n; ++i) {
    numerator += y[i] * (1 - d[i]);
    denominator += mu[i] * (1 - d[i]);
  }
  return numerator > denominator ? numerator / denominator : 1.0;
}

//' Estimate the indicators of excess zeros for a ZIP distribution.
//' 
//' Given counts and (estimated) Poisson expected value parameters and excess 
//' zero probabilities, this function estimates the indicator variable for an 
//' excess zero, for each count.
//' @param p A numeric vector, each element being the given/estimated 
//'    probability of an excess zero for the corresponding count in \code{y}.
//' @param mu A numeric vector, each element being the given/estimated Poisson 
//'    expected value parameter for the corresponding count in \code{y}. Of same 
//'    length as \code{p}.
//' @param y An integer vector containing the observed counts. Of same length as 
//'    \code{p}.
//' @return A numeric vector, of same length as the input vector \code{p}.
//' @keywords internal
// [[Rcpp::export]]
NumericVector estimate_d(const NumericVector& p, 
                         const NumericVector& mu, 
                         const NumericVector& y) {
  int n = y.size();
  NumericVector d(n);
  for (int i = 0; i < n; ++i) {
    if (y[i] > 0.0) {
      d[i] = 0.0;
    } else {
      d[i] = p[i] / (p[i] + (1 - p[i]) * exp(-mu[i]));
    }
  }
  return d;
}

//' Estimates the ZIP relative risk and excess zero indicators for a window.
//' 
//' For a single spatial or space-time window, this function uses the EM 
//' algorithm to estimate the relative risk and the excess zero indicators for 
//' counts assumed to be generated from a zero-inflated Poisson distribution.
//' @param p A numeric vector of the given/estimated excess zero probabilities 
//'    corresponding to each count.
//' @param mu A numeric vector of the given/estimated Poisson expected value
//'    parameters corresponding to each count. Of same length as \code{p}.
//' @param y An integer vector of the observed counts, of same length as 
//'    \code{p}.
//' @param tol A scalar between 0 and 1. It is the absolute tolerance criterion
//'    for the estimate of the excess zero indicator; convergence is reached when
//'    two successive elements in the sequence of estimates have an absolute 
//'    difference less than \code{tol}.
//' @return A list with two elements:
//' \describe{
//'   \item{q}{Scalar estimate of the relative risk.}
//'   \item{dstar}{Estimates of the excess zero indicator variables.}
//' }
//' @keywords internal
// [[Rcpp::export]]
List zip_em_estimates(const NumericVector& p,
                      const NumericVector& mu,
                      const NumericVector& y,
                      double tol = 0.01) {
  int n = y.size();
  NumericVector d(n);
  NumericVector d_prev = estimate_d(p, mu, y);
  double q = estimate_zip_relrisk(d_prev, mu, y);
  d = estimate_d(p, q * mu, y);
  double maxdiff = max(abs(d - d_prev));
  
  while (maxdiff > tol) {
    d_prev = d;
    q = estimate_zip_relrisk(d_prev, mu, y);
    d = estimate_d(p, q * mu, y);
    maxdiff = max(abs(d - d_prev));
  }
  return List::create(Rcpp::Named("q") = estimate_zip_relrisk(d, mu, y),
                      Rcpp::Named("dstar") = d);
}

//' Calculate a term in the sum of the logarithm of the ZIP window statistic.
//' 
//' This function calculates a term which appears in the sum of the logarithm of
//' the zero-inflated Poisson statistic for a given space-time window.
//' @param q Scalar; the relative risk.
//' @param p Numeric vector of excess zero probabilities.
//' @param dstar Numeric vector of estimates of the excess zero indicators, under 
//'    the alternative hypothesis of an outbreak. Of same length as \code{p}.
//' @param ddagger Numeric vector of estimates of the excess zero indicators, 
//'    under the null hypothesis of no outbreak. Of same length as \code{p}.
//' @param mu Numeric vector of given/estimated Poisson expected value 
//'    parameters. Of same length as \code{p}.
//' @param y Integer vector of observed counts. Of same length as \code{p}.
//' @return A numeric vector of same length as input vector \code{p}.
//' @keywords internal
// [[Rcpp::export]]
NumericVector zip_statistic_term(double q,
                                 const NumericVector& p,
                                 const NumericVector& dstar,
                                 const NumericVector& ddagger,
                                 const NumericVector& mu,
                                 const NumericVector& y) {
  return zip_statistic_logfactor(p, dstar, q * mu, y) - 
    zip_statistic_logfactor(p, ddagger, mu, y);
}

//' Calculate the ZIP statistic for a single space-time window.
//' 
//' Calculate the single-window statistic for the zero-inflated Poisson 
//' distribution using the EM algorithm.
//' @inheritParams zip_em_estimates
//' @param ... Named parameters passed to \code{\link{zip_em_estimates}}.
//' @return A scalar, the (logarithm of the) ZIP statistic.
//' @keywords internal
// [[Rcpp::export]]
double window_zip_statistic(const NumericVector& p,
                            const NumericVector& mu,
                            const NumericVector& y,
                            double tol = 0.01) {
  
  List em = zip_em_estimates(p, mu, y, tol);
  return sum(zip_statistic_term(em["q"], 
                                p, 
                                em["dstar"],
                                estimate_d(p, mu, y),
                                mu, 
                                y));
}


//' Calculate the ZIP window statistic over all durations, for a given zone.
//' 
//' This function calculates the zero-inflated Poisson statistic for a given 
//' spatial zone, for all durations considered.
//' @param duration An integer vector.
//' @inheritParams zip_em_estimates
//' @return A list with two elements:
//' \describe{
//'   \item{duration}{Vector of integers from 1 to \code{maxdur}.}
//'   \item{statistic}{Numeric vector containing the ZIP statistics corresponding
//'   to each duration, for the given spatial zone.}
//' }
//' @keywords internal
// [[Rcpp::export]]
List calc_zipstat_over_duration(const IntegerVector& duration,
                                const NumericVector& p,
                                const NumericVector& mu,
                                const NumericVector& y,
                                int maxdur,
                                double tol = 0.01) {
  IntegerVector dur(maxdur);
  NumericVector stat(maxdur);
  for (int i = 0; i < maxdur; ++i) {
    stat[i] = i + 1;
    stat[i] = window_zip_statistic(p[duration < i + 2], 
                                   mu[duration < i + 2], 
                                   y[duration < i + 2], 
                                   tol);
  }
  return List::create(Rcpp::Named("duration") = dur,
                      Rcpp::Named("statistic") = stat);
}
