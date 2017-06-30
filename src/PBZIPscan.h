#ifndef PBZIPSCAN_H
#define PBZIPSCAN_H

#include "USTscan.h"
#include "ZIPutility.h"
#include "scan_utility.h"

class PBZIPscan : public USTscan<arma::umat> {

public:
  PBZIPscan(const arma::umat& counts,
            const arma::mat& pop,
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
  arma::mat m_pop;
  double m_rel_tol;

  // Components of returned list
  arma::uvec m_zone_numbers;
  arma::uvec m_durations;
  arma::vec  m_scores;
  arma::vec  m_relrisk_in;
  arma::vec  m_relrisk_out;
  arma::vec  m_p;

  // Functions
  using store_ptr = void (PBZIPscan::*)(int storage_index, double score,
                                        double q_in, double q_out, double p,
                                        int zone_nr, int duration);
  store_ptr store;
  void store_max(int storage_index, double score, double q_in, double q_out,
                 double p, int zone_nr, int duration);
  void store_all(int storage_index, double score, double q_in, double q_out,
                 double p, int zone_nr, int duration);

  double pb_zip_relrisk(const int y_sum, const arma::vec& pop,
                        const arma::vec& d);

};

// Implementations -------------------------------------------------------------

inline PBZIPscan::PBZIPscan(const arma::umat& counts,
                            const arma::mat& pop,
                            const arma::uvec& zones,
                            const arma::uvec& zone_lengths,
                            const int num_locs,
                            const int num_zones,
                            const int max_dur,
                            const double rel_tol,
                            const bool store_everything)
  : USTscan(counts, zones, zone_lengths, num_locs, num_zones, max_dur),
    m_pop(pop),
    m_rel_tol(rel_tol) {

  int out_length;
  if (store_everything) {
    out_length = m_num_zones * m_max_dur;
    store = &PBZIPscan::store_all;
  } else {
    out_length = 1;
    store = &PBZIPscan::store_max;
  }

  m_zone_numbers.set_size(out_length);
  m_durations.set_size(out_length);
  m_scores.set_size(out_length);
  m_relrisk_in.set_size(out_length);
  m_relrisk_out.set_size(out_length);
  m_p.set_size(out_length);

  if (!store_everything) {
    m_scores[0] = -1.0;
  }
}

inline void PBZIPscan::calculate(const int storage_index,
                                 const int zone_nr,
                                 const int duration,
                                 const arma::uvec& current_zone,
                                 const arma::uvec& current_rows) {
  if (current_rows.n_elem == m_max_dur &&
      current_zone.n_elem == m_num_locs) {
    (this->*store)(storage_index, -1.0, -1.0, -1.0, -1.0, zone_nr + 1,
                   duration + 1);
  }


  Rcpp::Rcout << "storage_idx = "<< storage_index << std::endl;
  Rcpp::Rcout << "zone_nr = "<< zone_nr << std::endl;
  Rcpp::Rcout << "duration = "<< duration + 1 << std::endl;

  // Index vectors for outside current window
  arma::uvec out_zone = zone_complement(current_zone);
  arma::uvec out_rows = duration_complement(duration);
  Rcpp::Rcout << "out_zone = "<< out_zone.t() << std::endl;
  Rcpp::Rcout << "out_rows = "<< out_rows.t() << std::endl;

  // Extract counts and parameters as vectors
  arma::uvec y_in = arma::vectorise(m_counts.submat(current_rows,
                                                    current_zone));
  arma::vec  pop_in = arma::vectorise(m_pop.submat(current_rows, current_zone));

  arma::uvec y_out;
  arma::vec pop_out;

  // If current window has maximum duration
  if (duration + 1 == m_max_dur) {
    y_out = arma::vectorise(m_counts.cols(out_zone));
    pop_out = arma::vectorise(m_pop.cols(out_zone));
  // If current window covers all locations
  } else if (current_zone.n_elem == m_num_locs) {
    y_out = arma::vectorise(m_counts.rows(out_rows));
    pop_out = arma::vectorise(m_pop.rows(out_rows));
  } else {
    arma::uvec y_out = arma::vectorise(m_counts.submat(out_rows, out_zone));
    arma::vec  pop_out = arma::vectorise(m_pop.submat(out_rows, out_zone));
  }

  int n_in = y_in.n_elem;
  int n_out = y_out.n_elem;

  int y_sum_in = arma::accu(y_in);
  int y_sum_out = arma::accu(y_out);

  std::vector<int> zero_idx_in  = get_zero_indices(y_in);
  std::vector<int> zero_idx_out = get_zero_indices(y_out);

  // Initialization (Expectation step)
  arma::vec d_in(n_in);   d_in.fill(0.5);
  arma::vec d_out(n_out); d_out.fill(0.5);

  double p  = (arma::accu(d_in) + arma::accu(d_out))  / (n_in + n_out);

  // Run maximization step to get initial loglihood
  double q_in = pb_zip_relrisk(y_sum_in, pop_in, d_in);
  double q_out = pb_zip_relrisk(y_sum_out, pop_out, d_out);

  double loglik_old = zip_loglihood(y_in, pop_in, p * arma::ones(n_in), q_in)
                    + zip_loglihood(y_out, pop_out, p * arma::ones(n_out),
                                    q_out);
  double loglik_new;

  // Run EM algorithm
  double diff;
  int n_iterations = 1;

  do {
    // Expectation step
    for (const int& i : zero_idx_in) {
      d_in[i] = zip_zeroindic(pop_in[i], p, q_in);
    }

    for (const int& i : zero_idx_out) {
      d_out[i] = zip_zeroindic(pop_out[i], p, q_out);
    }

    p  = (arma::accu(d_in) + arma::accu(d_out))  / (n_in + n_out);


    // Maximization step
    q_in = pb_zip_relrisk(y_sum_in, pop_in, d_in);
    q_out = pb_zip_relrisk(y_sum_out, pop_out, d_out);

    // Update likelihood
    loglik_new = zip_loglihood(y_in, pop_in, p * arma::ones(n_in), q_in)
               + zip_loglihood(y_out, pop_out, p * arma::ones(n_out), q_out);
    diff = std::abs(exp(loglik_new - loglik_old) - 1.0);
    loglik_old = loglik_new;
    ++n_iterations;

  } while (diff > m_rel_tol);

  (this->*store)(storage_index, loglik_new, q_in, q_out, p, zone_nr + 1,
                 duration + 1);
}

inline Rcpp::DataFrame PBZIPscan::get_results() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")     = m_zone_numbers,
    Rcpp::Named("duration") = m_durations,
    Rcpp::Named("score")    = m_scores,
    Rcpp::Named("relrisk_in")  = m_relrisk_in,
    Rcpp::Named("relrisk_out")  = m_relrisk_out,
    Rcpp::Named("p")  = m_p);
}

inline void PBZIPscan::store_all(int storage_index, double score, double q_in,
                                 double p, double q_out, int zone_nr,
                                 int duration) {
  m_scores[storage_index]       = score;
  m_relrisk_in[storage_index]   = q_in;
  m_relrisk_out[storage_index]  = q_out;
  m_zone_numbers[storage_index] = zone_nr;
  m_durations[storage_index]    = duration;
  m_p[storage_index]            = p;
}

inline void PBZIPscan::store_max(int storage_index, double score, double q_in,
                                 double p, double q_out, int zone_nr,
                                 int duration) {
  if (score > m_scores[0]) {
    m_scores[0]       = score;
    m_relrisk_in[0]   = q_in;
    m_relrisk_out[0]  = q_out;
    m_zone_numbers[0] = zone_nr;
    m_durations[0]    = duration;
    m_p[0]            = p;
  }
}

inline double PBZIPscan::pb_zip_relrisk(const int y_sum, const arma::vec& pop,
                                        const arma::vec& d) {
  double denominator = 0.0;
  for (int i = 0; i < pop.n_elem; ++i) {
    denominator += pop[i] * (1.0 - d[i]);
  }
  return y_sum / denominator;
}


#endif
