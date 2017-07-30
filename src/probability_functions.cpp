#include <cmath>
#include "probability_functions.h"
using namespace Rcpp;


// Poisson distribution --------------------------------------------------------

// The log probability mass function of the Poisson distribution.
//
// The log probability mass function of the Poisson distribution, evaluated at
// \code{y}.
// @param y An integer.
// @param mu A positive scalar.
// @return A log-probability.
// @keywords internal
double poisson_lpmf(const double y, const double mu) {
  return y * log(mu) - lgamma(y + 1.0) - mu;
}

// Calculate the Poisson log-likelihood.
//
// Calculate the log-likelihood for the Poisson distribution.
// @param y A non-negative integer vector/matrix; the observed counts.
// @param mu A vector/matrix of positive scalars; the expected values of the
//   counts.
// @param q A scalar greater than or equal to 1; the relative risk.
// @return A non-positive scalar; the log-likelihood.
// @keywords internal
double poisson_loglihood(const arma::uvec &y,
                         const arma::vec &mu,
                         const double q) {
  double loglihood = 0.0;
  for (int i = 0; i < y.n_elem; ++i) {
    loglihood += poisson_lpmf(y(i), q * mu(i));
  }
  return loglihood;
}

// Negative binomial distribution ----------------------------------------------

// Draw an observation from the negative binomial distribution parametrized by
// mean \eqn{\mu} and overdispersion \eqn{\omega = 1 + \mu / \theta}.
int rnbinom2(const double mu, const double omega) {
  return (omega - 1.0 < 1e-9 ?
          R::rpois(mu) :
          R::rnbinom(mu / (omega - 1.0), 1.0 / omega));
}

// Zero-inflated Poisson distribution ------------------------------------------

// The log probability mass function of the ZIP distribution.
//
// The log probability mass function of the Poisson distribution, evaluated at
// \code{x}.
// @param y A non-negative integer; the observed count.
// @param mu A positive scalar; the expected value of the count.
// @param p A scalar between 0 and 1; the structural zero probability.
// @return A non-positive scalar; the loglihood contribution of the
//    observation.
// @keywords internal
double zip_lpmf(const int y, const double mu, const double p) {
  if (y == 0) {
    return log(p + (1 - p) * exp(-mu));
  } else {
    return log(1 - p) + y * log(mu) - lgamma(y + 1.0) - mu;
  }
}

// Calculate the ZIP log-likelihood.
//
// Calculate the (incomplete information) log-likelihood for the zero-inflated
// Poisson distribution.
// @param y A vector/matrix of non-negative integers; the observed counts.
// @param mu A vector/matrix of positive scalars; the expected values of the
//    counts.
// @param p A vector of scalars between 0 and 1; the structural zero
//    probabilities.
// @param q A scalar greater than or equal to 1; the relative risk.
// @return A non-positive scalar; the incomplete information ZIP loglihood.
// @keywords internal
double zip_loglihood(const arma::uvec &y,
                     const arma::vec &mu,
                     const arma::vec &p,
                     const double q) {
  double loglihood = 0.0;
  for (int i = 0; i < y.n_elem; ++i) {
    loglihood += zip_lpmf(y(i), q * mu(i), p(i));
  }
  return loglihood;
}

// Draw an observation from the ZIP distribution.
//
// Draw a sample (one observation) from the zero-inflated Poisson distribution.
// @param mu Scalar; The Poisson mean parameter.
// @param p Scalar; the structural zero probability.
// @return An integer.
// @keywords internal
int rzip(const double mu, const double p) {
  return (R::runif(0.0, 1.0) < p ? 0 : static_cast<int>(R::rpois(mu)));
}

