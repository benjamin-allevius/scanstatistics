#ifndef EBZIPSCAN_H
#define EBZIPSCAN_H

#include "USTscan.h"
#include "ZIPutility.h"
#include "scan_utility.h"

class EBZIPscan : public USTscan<arma::umat> {

public:
  EBZIPscan(const arma::umat& counts,
            const arma::mat& baselines,
            const arma::mat& probs,
            const arma::uvec& zones,
            const arma::uvec& zone_lengths,
            const int num_locs,
            const int num_zones,
            const int max_dur,
            const double rel_tol,
            const bool store_everything);
  void calculate(const int storage_index,
                 const int zone_nr,
                 const int duration,
                 const arma::uvec& current_zone,
                 const arma::uvec& current_rows);
  Rcpp::DataFrame get_results();

private:
  arma::mat m_baselines;
  arma::mat m_probs;
  double    m_rel_tol;

  // Components of returned list
  arma::uvec m_zone_numbers;
  arma::uvec m_durations;
  arma::vec  m_scores;
  arma::vec  m_relrisks;
  arma::uvec m_iterations;

  // Functions
  using store_ptr = void (EBZIPscan::*)(int storage_index, double score,
                                        double q, int n_iterations,
                                        int zone_nr, int duration);
  store_ptr store;
  void store_max(int storage_index, double score, double q, int n_iterations,
                 int zone_nr, int duration);
  void store_all(int storage_index, double score, double q, int n_iterations,
                 int zone_nr, int duration);
  double eb_zip_relrisk(const int y_sum, const arma::vec& mu,
                        const arma::vec& d);

};

// Implementations -------------------------------------------------------------

inline EBZIPscan::EBZIPscan(const arma::umat& counts,
                            const arma::mat& baselines,
                            const arma::mat& probs,
                            const arma::uvec& zones,
                            const arma::uvec& zone_lengths,
                            const int num_locs,
                            const int num_zones,
                            const int max_dur,
                            const double rel_tol,
                            const bool store_everything)
  : USTscan(counts, zones, zone_lengths, num_locs, num_zones, max_dur),
    m_baselines(baselines),
    m_probs(probs),
    m_rel_tol(rel_tol) {

  int out_length;
  if (store_everything) {
    out_length = m_num_zones * m_max_dur;
    store = &EBZIPscan::store_all;
  } else {
    out_length = 1;
    store = &EBZIPscan::store_max;
  }

  m_zone_numbers.set_size(out_length);
  m_durations.set_size(out_length);
  m_scores.set_size(out_length);
  m_relrisks.set_size(out_length);
  m_iterations.set_size(out_length);

  if (!store_everything) {
    m_scores[0] = -1.0;
  }
}

inline void EBZIPscan::calculate(const int storage_index,
                                 const int zone_nr,
                                 const int duration,
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
  std::vector<int> zero_idx = get_zero_indices(y);
  int y_sum = arma::accu(y);

  // Run EM algorithm
  double diff = 2.0 * m_rel_tol;
  int n_iterations = 0;
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

inline Rcpp::DataFrame EBZIPscan::get_results() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")     = m_zone_numbers,
    Rcpp::Named("duration") = m_durations,
    Rcpp::Named("score")    = m_scores,
    Rcpp::Named("relrisk")  = m_relrisks,
    Rcpp::Named("n_iter")   = m_iterations);
}

inline void EBZIPscan::store_all(int storage_index, double score, double q,
                                 int n_iterations, int zone_nr, int duration) {
  m_scores[storage_index]       = score;
  m_relrisks[storage_index]     = q;
  m_iterations[storage_index]   = n_iterations;
  m_zone_numbers[storage_index] = zone_nr;
  m_durations[storage_index]    = duration;
}

inline void EBZIPscan::store_max(int storage_index, double score, double q,
                                 int n_iterations, int zone_nr, int duration) {
  if (score > m_scores[0]) {
    m_scores[0]       = score;
    m_relrisks[0]     = q;
    m_iterations[0]   = n_iterations;
    m_zone_numbers[0] = zone_nr;
    m_durations[0]    = duration;
  }
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
inline double EBZIPscan::eb_zip_relrisk(const int y_sum, const arma::vec& mu,
                             const arma::vec& d) {
  double denominator = 0.0;
  for (int i = 0; i < mu.n_elem; ++i) {
    denominator += mu[i] * (1.0 - d[i]);
  }
  return std::max(1.0, y_sum / denominator);
}

#endif
