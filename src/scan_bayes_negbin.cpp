#include "scan_bayes_negbin.h"
#include "BGPscan.h"

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
                                 const arma::vec& m_probs) {
  BGPscan ob {counts, baselines, zones, zone_lengths, outbreak_prob,
              alpha_null, beta_null, alpha_alt, beta_alt, m_values, m_probs};
  ob.run_over_inc();
  return ob.get_results();
}

