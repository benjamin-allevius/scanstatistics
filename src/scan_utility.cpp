#include "scan_utility.h"

std::vector<int> get_zero_indices(const arma::uvec& v) {
  std::vector<int> zero_idx;
  for (int i = 0; i < v.n_elem; ++i) {
    if (v[i] == 0) zero_idx.push_back(i);
  }
  return zero_idx;
}

Rcpp::NumericVector armaVec2rcppVec(const arma::vec& v) {
  Rcpp::NumericVector x(v.size());
  for (int i = 0; i < v.size(); ++i) x.at(i) = v.at(i);
  return x;
}

arma::uvec rcppIVec2armaVec(const Rcpp::IntegerVector& v) {
  arma::uvec x(v.size());
  for (int i = 0; i < v.size(); ++i) x.at(i) = v.at(i);
  return x;
}
