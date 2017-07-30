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
std::vector<int> get_zero_indices(const arma::uvec& v);

// Convert an Rcpp NumericVector to an Armadillo vec.
arma::vec RcppVec_to_armaVec(const Rcpp::NumericVector& v);

// Convert an Armadillo vec to an Rcpp NumericVector
Rcpp::NumericVector armaVec2rcppVec(const arma::vec& v);

// Convert an Rcpp IntegerVector to an Armadillo uvec.
arma::uvec rcppIVec2armaVec(const Rcpp::IntegerVector& v);

// Comment: could not make template version of above 2 functions to work

#endif
