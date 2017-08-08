#ifndef SCAN_EB_ZIP_H
#define SCAN_EB_ZIP_H

#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]

//' Calculate the highest-value EB ZIP loglihood ratio statistic.
//'
//' Calculate the expectation-based ZIP loglihood ratio statistic for each zone
//' and duration, but only keep the zone and duration with the highest value
//' (the MLC). The estimate of the relative risk is also calculated, along with
//' the number of iterations the EM algorithm performed.
//' @param counts matrix (most recent timepoint in first row)
//' @param baselines matrix (most recent timepoint in first row)
//' @param probs matrix (most recent timepoint in first row)
//' @param zones integer vector (all zones concatenated; locations indexed from
//'    0 and up)
//' @param zone_lengths integer vector
//' @param rel_tol double
//' @param store_everything boolean
//' @param num_mcsim int
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
Rcpp::List scan_eb_zip_cpp(const arma::umat& counts,
                           const arma::mat& baselines,
                           const arma::mat& probs,
                           const arma::uvec& zones,
                           const arma::uvec& zone_lengths,
                           const double rel_tol,
                           const bool store_everything,
                           const arma::uword num_mcsim);


#endif
