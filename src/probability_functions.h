#ifndef PROBABILITY_FUNCTIONS_H
#define PROBABILITY_FUNCTIONS_H

#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]

// Poisson distribution
double poisson_lpmf(const double x, const double mu);
double poisson_loglihood(const arma::uvec &y,
                         const arma::vec &mu,
                         const double q = 1.0);

// Negative binomial distribution
int rnbinom2(const double mu, const double omega);

// ZIP distribution
double zip_lpmf(const int y, const double mu, const double p);
double zip_loglihood(const arma::uvec &y,
                     const arma::vec &mu,
                     const arma::vec &p,
                     const double q = 1.0);
int rzip(const double mu, const double p);


#endif
