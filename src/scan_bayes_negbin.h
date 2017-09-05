#ifndef SCAN_BAYES_NEGBIN_H
#define SCAN_BAYES_NEGBIN_H

#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]

//' Calculate the "Bayesian Spatial Scan Statistic" by Neill et al. (2006).
//' 
//' Calculate the "Bayesian Spatial Scan Statistic" by Neill et al. (2006),
//' adapted to a spatio-temporal setting. The scan statistic assumes that,
//' given the relative risk, the data follows a Poisson distribution. The 
//' relative risk is in turn assigned a Gamma distribution prior, yielding a 
//' negative binomial marginal distribution for the counts.
//' @param counts An integer matrix (most recent timepoint in first row).
//' @param baselines A matrix with positive entries (most recent timepoint in 
//'    first row).
//' @param zones An integer vector (all zones concatenated; locations indexed 
//'    from 0 and up).
//' @param zone_lengths An integer vector.
//' @param outbreak_prob A scalar; the probability of an outbreak (at any time,
//'    any place).
//' @param alpha_null A scalar; the shape parameter for the gamma distribution
//'    under the null hypothesis of no anomaly.
//' @param beta_null A scalar; the scale parameter for the gamma distribution
//'    under the null hypothesis of no anomaly.
//' @param alpha_alt A scalar; the shape parameter for the gamma distribution
//'    under the alternative hypothesis of an anomaly.
//' @param beta_alt A scalar; the scale parameter for the gamma distribution
//'    under the alternative hypothesis of an anomaly.
//' @param m_values A vector of possible values for the increase in the mean
//'    (and variance) of an anomalous count.
//' @param m_probs A vector of the prior probabilities of each value in 
//'    \code{m_values}.
//' @return A list with elements \code{observed} and \code{simulated}, each 
//'    being a data frame with columns:
//'    \describe{
//'      \item{zone}{The top-scoring zone (spatial component of MLC).}
//'      \item{duration}{The corresponding duration (time-length of MLC).}
//'      \item{score}{The value of the loglihood ratio statistic (the scan
//'                   statistic).}
//'    }
//' @export
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List scan_bayes_negbin_cpp(const arma::umat& counts,
                                 const arma::mat& baselines,
                                 const arma::uvec& zones,
                                 const arma::uvec& zone_lengths,
                                 const double outbreak_prob,
                                 const double alpha_null,
                                 const double beta_null,
                                 const double alpha_alt,
                                 const double beta_alt,
                                 const arma::vec& m_values,
                                 const arma::vec& m_probs);


#endif
