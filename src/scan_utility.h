#ifndef SCAN_UTILITY_H
#define SCAN_UTILITY_H

#include <vector>
#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]


//' Get indices of zero elements in a vector.
//' @param v An integer vector.
//' @return A vector with the indices of elements equal to zero in \code{v}.
//'    Indices start at zero.
//' @keywords internal
//' @export
// [[Rcpp::export]]
std::vector<int> get_zero_indices(arma::uvec v);

#endif
