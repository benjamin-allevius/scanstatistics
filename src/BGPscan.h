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
          const double alt_prob,
          const double alpha_null,
          const double beta_null,
          const double alpha_alt,
          const double beta_alt,
          const arma::vec& inc_vals,
          const arma::vec& inc_prior);

  Rcpp::DataFrame get_scan() override;
  Rcpp::List get_results();
  void run_over_inc();

private:
  arma::mat m_clprobs; // Conditional Log Probabilities, P(H_1(W) | D, m)
  arma::vec m_lprobs; // P(H_1(W) | D)
  arma::mat m_baselines;
  arma::mat m_baselines_orig;
  arma::uword m_total_count;
  double m_total_baseline;

  // Probabilities
  double m_alt_lprior;
  arma::vec m_alt_clposterior; // conditional (on inc_vals) log posterior
  double m_alt_lposterior;
  
  double m_null_lprior;
  arma::vec m_null_clposterior; 
  double m_null_lposterior;
  
  double m_window_lprior;
  
  arma::vec m_data_clprob;
  double m_data_lprob;
  
  // Increase of mean in case of an anomaly
  arma::uword m_inc_idx;
  arma::vec m_inc_vals;
  arma::vec m_inc_lprior;
  arma::vec m_inc_lposterior;
  
  // Matrices for the total anomaly probability for each location-time 
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
  
  Rcpp::List get_priors();
  Rcpp::List get_posteriors();
  Rcpp::DataFrame get_inc_posterior();
};

// Implementations -------------------------------------------------------------

inline BGPscan::BGPscan(const arma::umat& counts,
                        const arma::mat& baselines,
                        const arma::uvec& zones,
                        const arma::uvec& zone_lengths,
                        const double alt_prob,
                        const double alpha_null,
                        const double beta_null,
                        const double alpha_alt,
                        const double beta_alt,
                        const arma::vec& inc_vals,
                        const arma::vec& inc_prior)
  : USTscanbase(counts, zones, zone_lengths, true),
    m_baselines_orig(baselines),
    m_alt_lprior(std::log(alt_prob)),
    m_null_lprior(std::log(1.0 - alt_prob)),
    m_inc_idx(0),
    m_inc_vals(inc_vals),
    m_inc_lprior(arma::log(inc_prior)),
    m_alpha_null(alpha_null),
    m_beta_null(beta_null),
    m_alpha_alt(alpha_alt),
    m_beta_alt(beta_alt) {

  m_total_count = arma::accu(counts);
  m_total_baseline = arma::accu(baselines);

  m_counts = arma::cumsum(counts);
  m_baselines = arma::cumsum(baselines);
  
  // Set sizes
  m_clprobs.set_size(m_out_length, inc_vals.n_elem);
  m_lprobs.set_size(m_out_length);
  m_data_clprob.set_size(inc_vals.n_elem);
  m_alt_clposterior.set_size(inc_vals.n_elem);
  m_null_clposterior.set_size(inc_vals.n_elem);

  m_window_lprior = m_alt_lprior - std::log(m_num_zones * m_max_dur);
  
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

  score = log_prob(C, B, m_inc_vals(m_inc_idx) * m_alpha_alt, m_beta_alt) 
        + log_prob(m_total_count - C, m_total_baseline - B,
                   m_alpha_null, m_beta_null);

  (this->*store)(storage_index, score, zone_nr + 1, duration + 1);
}

inline void BGPscan::post_process() {
  // Calculate the probability of no anomaly
  // P(H_0 | D, m) \propto P(D | H_0, m) P(H_0 | m)
  double null_lprob = log_prob(m_total_count, m_total_baseline,
                               m_alpha_null, m_beta_null)
                    + m_null_lprior;
  
  // Calculate the marginal probability of the data, P(D | m)
  m_data_clprob(m_inc_idx) = 
    log_sum_exp(null_lprob, 
                m_window_lprior + log_sum_exp(m_clprobs.col(m_inc_idx)));
  
  // Normalize the posterior probabilities
  // P(H_0 | D, m) = P(D | H_0, m) P(H_0 | m) / P(D | m)
  m_null_clposterior(m_inc_idx) = null_lprob - m_data_clprob(m_inc_idx);
  
  // P(H_1(W) | D, m) = P(D | H_1(W), m) P(H_1(W) | m) / P(D | m)
  for (arma::uword i = 0; i < m_scores.n_elem; ++i) {
    m_clprobs(i, m_inc_idx) = m_clprobs(i, m_inc_idx) 
                               + m_window_lprior 
                               - m_data_clprob(m_inc_idx);
  }
  
  // P(H_1 | D, m)
  m_alt_clposterior(m_inc_idx) = 
    std::log(1.0 - std::exp(m_null_clposterior(m_inc_idx)));
}

inline void BGPscan::aggregate_scores(const arma::uword storage_index,
                                      const arma::uword zone_nr,
                                      const arma::uword duration,
                                      const arma::uvec& current_zone,
                                      const arma::uvec& current_rows) {
  for (arma::uword loc : current_zone) {
    m_spacetime_posteriors.at(duration, loc) += m_scores(storage_index);
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

inline void BGPscan::run_over_inc() {
  // Run the scan for each value in m_inc_vals
  for (arma::uword i = 0; i < m_inc_vals.n_elem; ++i, ++m_inc_idx) {
    run_scan(true);
  }
  
  // Calculate the marginal probability of the data ----------------------------
  // P(D) = \sum_m P(D | m) P(m)
  m_data_lprob = log_sum_exp(m_data_clprob + m_inc_lprior);
  
  // P(m | D) = P(D | m) P(m) / P(D)
  m_inc_lposterior = m_data_clprob + m_inc_lprior - m_data_lprob;
  
  // Calculate the posterior probability for each space-time window ------------
  // P(H_1(W) | D) = \sum_m P(H_1(W) | D, m) P(m | D)
  for (arma::uword i = 0; i < m_scores.n_elem; ++i) {
    m_lprobs(i) = log_sum_exp(m_clprobs.row(i).t() + m_inc_lposterior);
  }
  m_scores = arma::exp(m_lprobs);
  
  // Calculate the null and alternative posterior probabilities ----------------
  m_null_lposterior = log_sum_exp(m_null_clposterior + m_inc_lposterior);
  
  m_alt_lposterior = std::log(1.0 - std::exp(m_null_lposterior));
  
  // Calculate posterior probabilities for each location and location-time
  score_locations();
  m_location_posteriors = (arma::sum(m_spacetime_posteriors)).t();
}


// Storage functions -----------------------------------------------------------

inline void BGPscan::store_all(arma::uword storage_index, double score,
                               arma::uword zone_nr, arma::uword duration) {
  m_clprobs(storage_index, m_inc_idx) = score;
  if (m_inc_idx == 0) {
    m_zone_numbers(storage_index) = zone_nr;
    m_durations(storage_index)    = duration;
  }
}

// Retrieval functions ---------------------------------------------------------

inline Rcpp::DataFrame BGPscan::get_scan() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")          = m_zone_numbers,
    Rcpp::Named("duration")      = m_durations,
    Rcpp::Named("log_posterior") = m_lprobs,
    Rcpp::Named("log_bayes_factor")  = m_lprobs - m_null_lposterior);
}

inline Rcpp::DataFrame BGPscan::get_inc_posterior() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("inc_value") = m_inc_vals,
    Rcpp::Named("inc_posterior") = arma::exp(m_inc_lposterior));
}

inline Rcpp::List BGPscan::get_priors() {
  return Rcpp::List::create(
    Rcpp::Named("null_prior")   = std::exp(m_null_lprior),
    Rcpp::Named("alt_prior")    = std::exp(m_alt_lprior),
    Rcpp::Named("inc_prior")    = arma::exp(m_inc_lprior).t(),
    Rcpp::Named("window_prior") = std::exp(m_window_lprior)
  );
}

inline Rcpp::List BGPscan::get_posteriors() {
  return Rcpp::List::create(
    Rcpp::Named("null_posterior")        = std::exp(m_null_lposterior),
    Rcpp::Named("alt_posterior")         = std::exp(m_alt_lposterior),
    Rcpp::Named("inc_posterior")         = get_inc_posterior(),
    Rcpp::Named("window_posteriors")     = get_scan(),
    Rcpp::Named("space_time_posteriors") = m_spacetime_posteriors,
    Rcpp::Named("location_posteriors")   = m_location_posteriors.t()
  );
}

inline Rcpp::List BGPscan::get_results() {
  return Rcpp::List::create(
    Rcpp::Named("priors")             = get_priors(),
    Rcpp::Named("posteriors")         = get_posteriors(),
    Rcpp::Named("marginal_data_prob") = exp(m_data_lprob)
  );
}

#endif
