#ifndef EBZIPSCAN_H
#define EBZIPSCAN_H

#include "USTscan.h"
#include "ZIPutility.h"
#include "scan_utility.h"

class EBZIPscan : public USTscan<arma::umat, arma::uword> {

public:
  EBZIPscan(const arma::umat& counts,
            const arma::mat& baselines,
            const arma::mat& probs,
            const arma::uvec& zones,
            const arma::uvec& zone_lengths,
            const double rel_tol,
            const bool store_everything,
            const arma::uword num_mcsim);

  Rcpp::DataFrame get_scan()  override;
  Rcpp::DataFrame get_mcsim() override;

private:
  arma::mat m_baselines;
  arma::mat m_probs;
  double    m_rel_tol;

  // Values calculated on observed data
  arma::vec  m_relrisks;
  arma::uvec m_iterations;

  // Values calculated on simulated data
  arma::vec  sim_relrisks;
  arma::uvec sim_iterations;

  // Functions
  void calculate(const arma::uword storage_index,
                 const arma::uword zone_nr,
                 const arma::uword duration,
                 const arma::uvec& current_zone,
                 const arma::uvec& current_rows) override;
  arma::uword draw_sample(arma::uword row, arma::uword col) override;
  void set_sim_store_fun() override;

  using store_ptr = void (EBZIPscan::*)(arma::uword storage_index, double score,
                                        double q, arma::uword n_iterations,
                                        arma::uword zone_nr, arma::uword duration);
  store_ptr store;
  void store_max(arma::uword storage_index, double score, double q, arma::uword n_iterations,
                 arma::uword zone_nr, arma::uword duration);
  void store_all(arma::uword storage_index, double score, double q, arma::uword n_iterations,
                 arma::uword zone_nr, arma::uword duration);
  void store_sim(arma::uword storage_index, double score, double q, arma::uword n_iterations,
                 arma::uword zone_nr, arma::uword duration);
  double eb_zip_relrisk(const arma::uword y_sum, const arma::vec& mu,
                        const arma::vec& d);

};

// Implementations -------------------------------------------------------------

inline EBZIPscan::EBZIPscan(const arma::umat& counts,
                            const arma::mat& baselines,
                            const arma::mat& probs,
                            const arma::uvec& zones,
                            const arma::uvec& zone_lengths,
                            const double rel_tol,
                            const bool store_everything,
                            const arma::uword num_mcsim)
  : USTscan(counts, zones, zone_lengths, store_everything, num_mcsim),
    m_baselines(baselines),
    m_probs(probs),
    m_rel_tol(rel_tol) {

  store = (m_store_everything ? &EBZIPscan::store_all : &EBZIPscan::store_max);

  // Reserve sizes for values calculated on observed data
  m_relrisks.set_size(m_out_length);
  m_iterations.set_size(m_out_length);

  // Reserve sizes for values calculated on simulated
  sim_relrisks.set_size(m_num_mcsim);
  sim_iterations.set_size(m_num_mcsim);
}

// Workhorse functions ---------------------------------------------------------

inline void EBZIPscan::calculate(const arma::uword storage_index,
                                 const arma::uword zone_nr,
                                 const arma::uword duration,
                                 const arma::uvec& current_zone,
                                 const arma::uvec& current_rows) {
  // Extract counts and parameters as vectors
  arma::uvec y = arma::vectorise(m_counts.submat(current_rows,
                                                 current_zone));
  arma::vec mu = arma::vectorise(m_baselines.submat(current_rows,
                                                    current_zone));
  arma::vec  p = arma::vectorise(m_probs.submat(current_rows,
                                                current_zone));

  arma::vec d_hat = arma::zeros(y.n_elem); // Structural zero estimates
  double q_hat = 1.0;                      // Relative risk estimate

  double loglik_null = zip_loglihood(y, mu, p, 1.0);
  double loglik_old = loglik_null;
  double loglik_new;

  // Store indices of zero counts; only loop through these when estimating
  // structural zero indicators. Also compute the sum of all counts.
  std::vector<arma::uword> zero_idx = get_zero_indices(y);
  arma::uword y_sum = arma::accu(y);

  // Run EM algorithm
  double diff = 2.0 * m_rel_tol;
  arma::uword n_iterations = 0;
  while (diff > m_rel_tol) {
    n_iterations += 1;

    // Expectation-step
    for (const int& i : zero_idx) {
      d_hat[i] = zip_zeroindic(mu[i], p[i], q_hat);
    }
    // Maximization-step
    q_hat = eb_zip_relrisk(y_sum, mu, d_hat);

    // Update likelihood
    loglik_new = zip_loglihood(y, mu, p, q_hat);
    diff = std::abs(exp(loglik_new - loglik_old) - 1.0);
    loglik_old = loglik_new;
  }

  (this->*store)(storage_index, loglik_new - loglik_null, q_hat, n_iterations,
                 zone_nr + 1, duration + 1);
}

// Estimate the relative risk for the ZIP distribution.
//
// Estimate the relative risk for the ZIP distribution.
// @param y_sum A non-negative integer; the sum of the observed counts.
// @param mu A vector of positive numbers; the expected values of the counts or
//    the corresponding population.
// @param d A vector of (estimates of) the structural zero indicators.
// @return A scalar; the relative risk.
// @keywords internal
inline double EBZIPscan::eb_zip_relrisk(const arma::uword y_sum, const arma::vec& mu,
                             const arma::vec& d) {
  double denominator = 0.0;
  for (arma::uword i = 0; i < mu.n_elem; ++i) {
    denominator += mu[i] * (1.0 - d[i]);
  }
  return std::max(1.0, y_sum / denominator);
}

inline arma::uword EBZIPscan::draw_sample(arma::uword row, arma::uword col) {
  return rzip(m_baselines.at(row, col), m_probs.at(row, col));
}

// Storage functions -----------------------------------------------------------

inline void EBZIPscan::store_all(arma::uword storage_index, double score, double q,
                                 arma::uword n_iterations, arma::uword zone_nr, arma::uword duration) {
  m_scores[storage_index]       = score;
  m_relrisks[storage_index]     = q;
  m_iterations[storage_index]   = n_iterations;
  m_zone_numbers[storage_index] = zone_nr;
  m_durations[storage_index]    = duration;
}

inline void EBZIPscan::store_max(arma::uword storage_index, double score, double q,
                                 arma::uword n_iterations, arma::uword zone_nr, arma::uword duration) {
  if (score > m_scores[0]) {
    m_scores[0]       = score;
    m_relrisks[0]     = q;
    m_iterations[0]   = n_iterations;
    m_zone_numbers[0] = zone_nr;
    m_durations[0]    = duration;
  }
}

inline void EBZIPscan::store_sim(arma::uword storage_index, double score, double q,
                                 arma::uword n_iterations, arma::uword zone_nr, arma::uword duration) {
  if (score > sim_scores[m_mcsim_index]) {
    sim_scores[m_mcsim_index]       = score;
    sim_relrisks[m_mcsim_index]     = q;
    sim_iterations[m_mcsim_index]   = n_iterations;
    sim_zone_numbers[m_mcsim_index] = zone_nr;
    sim_durations[m_mcsim_index]    = duration;
  }
}

inline void EBZIPscan::set_sim_store_fun() {
  store = &EBZIPscan::store_sim;
}

// Retrieval functions ---------------------------------------------------------

inline Rcpp::DataFrame EBZIPscan::get_scan() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")     = m_zone_numbers,
    Rcpp::Named("duration") = m_durations,
    Rcpp::Named("score")    = m_scores,
    Rcpp::Named("relrisk")  = m_relrisks,
    Rcpp::Named("n_iter")   = m_iterations);
}

inline Rcpp::DataFrame EBZIPscan::get_mcsim() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")     = sim_zone_numbers,
    Rcpp::Named("duration") = sim_durations,
    Rcpp::Named("score")    = sim_scores,
    Rcpp::Named("relrisk")  = sim_relrisks,
    Rcpp::Named("n_iter")   = sim_iterations);
}

#endif
