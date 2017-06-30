#ifndef ZIPUTILITY_H
#define ZIPUTILITY_H

#include <cmath>
#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]


// Estimate a structural zero indicator for the ZIP distribution.
//
// Estimate a structural zero indicator for the ZIP distribution.
// @param mu A positive scalar; the expected values of the counts or
//    the corresponding population.
// @param p A positive scalar; the structural zero probability.
// @param q A positive scalar; the relative risk.
// @return A scalar; the estimate of the structural zero indicator.
// @keywords internal
inline double zip_zeroindic(const double mu, const double p, const double q) {
  return p / (p + (1 - p) * std::exp(-q * mu));
}

#endif
