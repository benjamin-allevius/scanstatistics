#ifndef SCAN_EB_POI_H
#define SCAN_EB_POI_H

#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]

//' Calculate the expecation-based Poisson scan statistic.
//' 
//' Calculate the expectation-based Poisson scan statistic and Monte Carlo 
//' replicates.
//' @param counts An integer matrix (most recent timepoint in first row).
//' @param baselines A matrix with positive entries (most recent timepoint in 
//'    first row).
//' @param zones An integer vector (all zones concatenated; locations indexed 
//'    from 0 and up).
//' @param zone_lengths An integer vector.
//' @param store_everything A boolean.
//' @param num_mcsim An integer.
//' @return A list with elements \code{observed} and \code{simulated}, each 
//'    being a data frame with columns:
//'    \describe{
//'      \item{zone}{The top-scoring zone (spatial component of MLC).}
//'      \item{duration}{The corresponding duration (time-length of MLC).}
//'      \item{score}{The value of the loglihood ratio statistic (the scan
//'                   statistic).}
//'      \item{relrisk_in}{The estimated relative risk inside.}
//'      \item{relrisk_in}{The estimated relative risk outside.}
//'    }
//' @export
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List scan_eb_poisson_cpp(const arma::umat& counts,
                               const arma::mat& baselines,
                               const arma::uvec& zones,
                               const arma::uvec& zone_lengths,
                               const bool store_everything,
                               const arma::uword num_mcsim);


#endif
