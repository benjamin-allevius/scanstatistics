#ifndef BGPSCAN_H
#define BGPSCAN_H

#include "USTscan.h"
#include "scan_utility.h"

class BGPscan : public USTscanbase<arma::umat, arma::uword> {

public:
  BGPscan(const arma::umat& counts,
          const arma::mat& baselines,
          const arma::uvec& zones,
          const arma::uvec& zone_lengths,
          const double outbreak_prob,
          const double alpha_null,
          const double beta_null,
          const double alpha_alt,
          const double beta_alt);

  Rcpp::DataFrame get_scan() override;
  Rcpp::List get_results();
  Rcpp::List get_log_results();
  
  friend Rcpp::List scan_bayes_negbin_cpp(const arma::umat& counts,
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

private:
  arma::vec m_log_scores;
  arma::mat m_baselines;
  arma::mat m_baselines_orig;
  arma::uword m_total_count;
  double m_total_baseline;

  // Probabilities
  double m_log_outbreak_prior;
  double m_log_outbreak_posterior;
  
  double m_log_null_prior;
  double m_log_null_posterior;
  
  double m_log_window_prior;
  double m_log_data_prob;
  
  // Matrices for the total outbreak probability for each location-time 
  // combination and each location
  arma::mat m_spacetime_posteriors;
  arma::vec m_location_posteriors;

  // Relative risk gamma distribution parameters
  double m_alpha_null;
  double m_beta_null;
  double m_alpha_alt;
  double m_beta_alt;

  // Functions
  void calculate(const arma::uword storage_index,
                 const arma::uword zone_nr,
                 const arma::uword duration,
                 const arma::uvec& current_zone,
                 const arma::uvec& current_rows) override;
  
  void aggregate_scores(const arma::uword storage_index,
                        const arma::uword zone_nr,
                        const arma::uword duration,
                        const arma::uvec& current_zone,
                        const arma::uvec& current_rows) override;

  using store_ptr = void (BGPscan::*)(arma::uword storage_index, double score,
                          arma::uword zone_nr, arma::uword duration);
  store_ptr store;
  void store_all(arma::uword storage_index, double score, arma::uword zone_nr,
                 arma::uword duration);

  void post_process() override;

  double log_prob(const arma::uword C, const double B, const double alpha,
                  const double beta);

};

// Implementations -------------------------------------------------------------

inline BGPscan::BGPscan(const arma::umat& counts,
                        const arma::mat& baselines,
                        const arma::uvec& zones,
                        const arma::uvec& zone_lengths,
                        const double outbreak_prob,
                        const double alpha_null,
                        const double beta_null,
                        const double alpha_alt,
                        const double beta_alt)
  : USTscanbase(counts, zones, zone_lengths, true),
    m_baselines_orig(baselines),
    m_log_outbreak_prior(std::log(outbreak_prob)),
    m_log_null_prior(std::log(1.0 - outbreak_prob)),
    m_alpha_null(alpha_null),
    m_beta_null(beta_null),
    m_alpha_alt(alpha_alt),
    m_beta_alt(beta_alt) {

  m_total_count = arma::accu(counts);
  m_total_baseline = arma::accu(baselines);

  m_counts = arma::cumsum(counts);
  m_baselines = arma::cumsum(baselines);
  
  m_log_scores.set_size(m_out_length);

  m_log_window_prior = m_log_outbreak_prior
                     - (std::log(m_num_zones) + std::log(m_max_dur));
  
  m_spacetime_posteriors.zeros(m_max_dur, m_num_locs);
  m_location_posteriors.zeros(m_num_locs);

  store = &BGPscan::store_all;

}

// Workhorse functions ---------------------------------------------------------

inline void BGPscan::calculate(const arma::uword storage_index,
                               const arma::uword zone_nr,
                               const arma::uword duration,
                               const arma::uvec& current_zone,
                               const arma::uvec& current_rows) {
  arma::uword C;
  double B, score;

  arma::uvec row_idx = current_rows.tail(1);

  // Counts and baselines are already aggregated
  C = arma::accu(m_counts.submat(row_idx, current_zone));
  B = arma::accu(m_baselines.submat(row_idx, current_zone));

  score = log_prob(C, B, m_alpha_alt, m_beta_alt) 
        + log_prob(m_total_count - C, m_total_baseline - B,
                   m_alpha_null, m_beta_null);

  (this->*store)(storage_index, score, zone_nr + 1, duration + 1);
}

inline void BGPscan::post_process() {
  // Calculate the probability of no outbreak, P(H_0 | D) = P(D | H_0) P(H_0)
  double null_log_prob = log_prob(m_total_count, m_total_baseline,
                                  m_alpha_null, m_beta_null)
                       + m_log_null_prior;
  
  // Calculate the marginal probability of the data, P(D)
  m_log_data_prob = log_sum_exp(null_log_prob, 
                                m_log_window_prior + log_sum_exp(m_log_scores));
  
  // Normalize the posterior probabilities
  m_log_null_posterior = null_log_prob - m_log_data_prob;
  
  for (arma::uword i = 0; i < m_scores.n_elem; ++i) {
    m_log_scores.at(i) = m_log_scores.at(i) 
                       + m_log_window_prior 
                       - m_log_data_prob;
  }
  m_scores = arma::exp(m_log_scores);
  
  m_log_outbreak_posterior = std::log(1.0 - std::exp(m_log_null_posterior));
  
  // 
  score_locations();
  m_location_posteriors = (arma::sum(m_spacetime_posteriors)).t();
}

inline void BGPscan::aggregate_scores(const arma::uword storage_index,
                                      const arma::uword zone_nr,
                                      const arma::uword duration,
                                      const arma::uvec& current_zone,
                                      const arma::uvec& current_rows) {
  for (arma::uword loc : current_zone) {
    m_spacetime_posteriors.at(duration, loc) += m_scores.at(storage_index);
  }
}


inline double BGPscan::log_prob(const arma::uword C, const double B,
                                const double alpha, const double beta) {
  return   alpha * std::log(beta) 
         + std::lgamma(alpha + C) 
         - (alpha + C) * std::log(beta + B) 
         - std::lgamma(alpha)
         - std::lgamma(1.0 + C);
}


// Storage functions -----------------------------------------------------------

inline void BGPscan::store_all(arma::uword storage_index, double score,
                               arma::uword zone_nr, arma::uword duration) {
  m_log_scores[storage_index]   = score;
  m_zone_numbers[storage_index] = zone_nr;
  m_durations[storage_index]    = duration;
}

// Retrieval functions ---------------------------------------------------------

inline Rcpp::DataFrame BGPscan::get_scan() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")      = m_zone_numbers,
    Rcpp::Named("duration")  = m_durations,
    Rcpp::Named("posterior") = m_scores);
}

inline Rcpp::List BGPscan::get_results() {
  return Rcpp::List::create(
    Rcpp::Named("null_prior")            = std::exp(m_log_null_prior),
    Rcpp::Named("null_posterior")        = std::exp(m_log_null_posterior),
    Rcpp::Named("outbreak_prior")        = std::exp(m_log_outbreak_prior),
    Rcpp::Named("outbreak_posterior")    = std::exp(m_log_outbreak_posterior),
    Rcpp::Named("window_prior")          = std::exp(m_log_window_prior),
    Rcpp::Named("window_posteriors")     = get_scan(),
    Rcpp::Named("space_time_posteriors") = m_spacetime_posteriors,
    Rcpp::Named("location_posteriors")   = m_location_posteriors,
    Rcpp::Named("marginal_data_prob")    = exp(m_log_data_prob)
  );
}

inline Rcpp::List BGPscan::get_log_results() {
  return Rcpp::List::create(
    Rcpp::Named("null_prior")            = m_log_null_prior,
    Rcpp::Named("null_posterior")        = m_log_null_posterior,
    Rcpp::Named("outbreak_prior")        = m_log_outbreak_prior,
    Rcpp::Named("outbreak_posterior")    = m_log_outbreak_posterior,
    Rcpp::Named("window_prior")          = m_log_window_prior,
    Rcpp::Named("space_time_posteriors") = m_spacetime_posteriors,
    Rcpp::Named("location_posteriors")   = m_location_posteriors,
    Rcpp::Named("marginal_data_prob")    = m_log_data_prob,
    Rcpp::Named("window_zone")           = m_zone_numbers,
    Rcpp::Named("window_duration")       = m_durations,
    Rcpp::Named("window_log_posterior")  = m_log_scores
  );
}

#endif
