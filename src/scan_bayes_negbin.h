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
//' @param inc_values A vector of possible values for the increase in the mean
//'    (and variance) of an anomalous count.
//' @param inc_probs A vector of the prior probabilities of each value in 
//'    \code{inc_values}.
//' @return A list with elements \code{priors} (list), \code{posteriors} (list), 
//'    and \code{marginal_data_prob} (scalar). The list \code{priors} has 
//'    elements
//'    \describe{
//'      \item{null_prior}{The prior probability of no anomaly.}
//'      \item{alt_prior}{The prior probability of an anomaly.}
//'      \item{inc_prior}{A vector (matrix with 1 row) of prior probabilities
//'                       of each value in the argument \code{m_values}.}
//'      \item{window_prior}{The prior probability of an outbreak in any of the
//'                          space-time windows.}
//'    }
//'    The list \code{posteriors} has elements
//'    \describe{
//'      \item{null_posterior}{The posterior probability of no anomaly.}
//'      \item{alt_posterior}{The posterior probability of an anomaly.}
//'      \item{inc_posterior}{A data frame with columns \code{inc_values} and
//'                           \code{inc_posterior}.}
//'      \item{window_posteriors}{A data frame with columns \code{zone}, 
//'                               \code{duration}, \code{log_posterior} and 
//'                               \code{log_bayes_factor}, each row 
//'                               corresponding to a space-time window.}
//'      \item{space_time_posteriors}{A matrix with the posterior anomaly 
//'                                   probability of each location-time 
//'                                   combination.}
//'      \item{location_posteriors}{A vector (matrix with 1 row) with the 
//'                                 posterior probability of an anomaly at each
//'                                 location.}
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
                                 const arma::vec& inc_values,
                                 const arma::vec& inc_probs);


#endif
