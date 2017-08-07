#include "scan_utility.h"

std::vector<int> get_zero_indices(const arma::uvec& v) {
  std::vector<int> zero_idx;
  for (int i = 0; i < v.n_elem; ++i) {
    if (v[i] == 0) zero_idx.push_back(i);
  }
  return zero_idx;
}


// Rcpp vector to armadillo vectors --------------------------------------------

arma::vec NumericVector2vec(const Rcpp::NumericVector& v) {
  arma::vec x(v.size());
  for (int i = 0; i < v.size(); ++i) x.at(i) = v.at(i);
  return x;
}

arma::uvec IntegerVector2uvec(const Rcpp::IntegerVector& v) {
  arma::uvec x(v.size());
  for (int i = 0; i < v.size(); ++i) x.at(i) = v.at(i);
  return x;
}

// armadillo vectors to Rcpp vectors -------------------------------------------

Rcpp::NumericVector vec2NumericVector(const arma::vec& v) {
  Rcpp::NumericVector x(v.size());
  for (int i = 0; i < v.size(); ++i) x.at(i) = v.at(i);
  return x;
}

Rcpp::IntegerVector uvec2IntegerVector(const arma::uvec& v) {
  Rcpp::IntegerVector x(v.size());
  for (int i = 0; i < v.size(); ++i) x.at(i) = v.at(i);
  return x;
}

//

arma::umat expand_matrix(const arma::umat& A) {
  arma::umat res(static_cast<arma::uword>(arma::accu(A)), 2);
  arma::uword index = 0;
  for (arma::uword j = 0; j < A.n_cols; ++j) {
    for (arma::uword i = 0; i < A.n_rows; ++i) {
      for (arma::uword k = 0; k < A(i, j); ++k) {
        res(index, 0) = i;
        res(index, 1) = j;
        ++index;
      }
    }
  }
  return res;
}

arma::umat contract_matrix(const arma::umat& A, 
                           arma::uword nr, arma::uword nc) {
  arma::umat out(nr, nc, arma::fill::zeros);
  for (arma::uword i = 0; i < A.n_rows; ++i) {
    ++out(A(i, 0), A(i, 1)); 
  }
  return out;
}

// Permute column 0 using Fisher-Yates algorithm
arma::uvec shuffle_time_counts(const arma::uvec& v) {
  arma::uvec res(v);
  arma::uword k = 0;
  for (arma::uword i = res.n_elem - 1; k < res.n_elem; --i, ++k) {
    arma::uword j = static_cast<arma::uword>(R::runif(0, static_cast<double>(i)));
    arma::uword tmp = res.at(i);
    res.at(i) = res.at(j);
    res.at(j) = tmp;
  }
  return res;
}

