#ifndef SCAN_UTILITY_H
#define SCAN_UTILITY_H

#include <vector>
#include <cmath>
#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]


//' Get indices of zero elements in a vector.
//' @param v An integer vector.
//' @return A vector with the indices of elements equal to zero in \code{v}.
//'    Indices start at zero.
//' @keywords internal
//' @export
// [[Rcpp::export]]
std::vector<arma::uword> get_zero_indices(const arma::uvec& v);

// Rcpp vector to armadillo vectors --------------------------------------------

arma::vec NumericVector2vec(const Rcpp::NumericVector& v);

arma::uvec IntegerVector2uvec(const Rcpp::IntegerVector& v);

// armadillo vectors to Rcpp vectors -------------------------------------------

// Convert an armadillo vec to an Rcpp NumericVector
Rcpp::NumericVector vec2NumericVector(const arma::vec& v);

// Convert an armadillo vec to an Rcpp NumericVector
Rcpp::IntegerVector uvec2IntegerVector(const arma::uvec& v);


// Comment: could not make template version of above 2 functions to work

// Functions for permuting the counts in an integer matrix while preserving the
// row and column marginal sums.
// integer matrix --> expand_matrix --> permute one column --> contract matrix
arma::umat expand_matrix(const arma::umat& A);

arma::umat contract_matrix(const arma::umat& A, arma::uword nr, arma::uword nc);

// Permute using Fisher-Yates algorithm
arma::uvec shuffle_time_counts(const arma::uvec& v);

//' Permute the entries of the matrix, preserving row and column marginals.
//' 
//' Permute the entries of the matrix, preserving row and column marginals.
//' @param A An integer matrix.
//' @return An integer matrix.
//' @keywords internal
// [[Rcpp::export]]
arma::umat permute_matrix(const arma::umat& A);

double log_sum_exp(const arma::vec& v, 
                   const double start_val, 
                   const double max_val);

double log_sum_exp(const arma::vec& v, const double max_val);

double log_sum_exp(const arma::vec& v);

double log_sum_exp(const double a, const double b);

#endif
