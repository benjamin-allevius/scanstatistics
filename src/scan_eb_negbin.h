#ifndef SCAN_EB_NEGBIN_H
#define SCAN_EB_NEGBIN_H

#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]

//' Calculate the expectation-based negative binomial scan statistic.
//'
//' Calculate the expectation-based negative binomial scan statistic and Monte
//' Carlo replicates.
//' @param counts Integer matrix (most recent timepoint in first row)
//' @param baselines Matrix (most recent timepoint in first row)
//' @param overdisp Matrix (most recent timepoint in first row)
//' @param zones Integer vector (all zones concatenated; locations indexed from
//'    0 and up)
//' @param zone_lengths Integer vector
//' @param store_everything Boolean
//' @param num_mcsim Integer
//' @param score_hotspot Boolean
//' @return A list with elements \code{observed} and \code{simulated}, each 
//'    being a data frame with columns:
//'    \describe{
//'      \item{zone}{The top-scoring zone (spatial component of MLC).}
//'      \item{duration}{The corresponding duration (time-length of MLC).}
//'      \item{score}{The value of the loglihood ratio statistic (the scan
//'                   statistic).}
//'      \item{relrisk}{The estimated relative risk.}
//'      \item{n_iter}{The number of iterations performed by the EM algorithm.}
//'    }
//' @export
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List scan_eb_negbin_cpp(const arma::umat& counts,
                              const arma::mat& baselines,
                              const arma::mat& overdisp,
                              const arma::uvec& zones,
                              const arma::uvec& zone_lengths,
                              const bool store_everything,
                              const arma::uword num_mcsim,
                              const bool score_hotspot);


#endif
