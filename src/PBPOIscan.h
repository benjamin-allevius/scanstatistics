#ifndef PBPOISCAN_H
#define PBPOISCAN_H

#include "USTscan.h"
#include "scan_utility.h"
#include <RcppArmadilloExtensions/rmultinom.h>

class PBPOIscan : public USTscan<arma::umat, int> {

public:
  PBPOIscan(const arma::umat& counts,
            const arma::mat& baselines,
            const arma::uvec& zones,
            const arma::uvec& zone_lengths,
            const bool store_everything,
            const int num_mcsim);

  Rcpp::DataFrame get_scan()  override;
  Rcpp::DataFrame get_mcsim() override;

private:
  arma::mat m_baselines;
  arma::mat m_baselines_orig; // used for simulation
  int m_total_count;

  // Values calculated on observed data
  arma::vec  m_relrisk_in;
  arma::vec  m_relrisk_out;

  // Values calculated on simulated data
  arma::vec  sim_relrisk_in;
  arma::vec  sim_relrisk_out;

  // Functions
  void calculate(const int storage_index,
                 const int zone_nr,
                 const int duration,
                 const arma::uvec& current_zone,
                 const arma::uvec& current_rows) override;
  void simulate_counts() override;
  void set_sim_store_fun() override;
  int draw_sample(arma::uword row, arma::uword col) override;

  using store_ptr = void (PBPOIscan::*)(int storage_index, double score,
                                        double q_in, double q_out, int zone_nr,
                                        int duration);
  store_ptr store;
  void store_max(int storage_index, double score, double q_in, double q_out,
                 int zone_nr, int duration);
  void store_all(int storage_index, double score, double q_in, double q_out,
                 int zone_nr, int duration);
  void store_sim(int storage_index, double score, double q_in, double q_out,
                 int zone_nr, int duration);

};

// Implementations -------------------------------------------------------------

inline PBPOIscan::PBPOIscan(const arma::umat& counts,
                            const arma::mat& baselines,
                            const arma::uvec& zones,
                            const arma::uvec& zone_lengths,
                            const bool store_everything,
                            const int num_mcsim)
  : USTscan(counts, zones, zone_lengths, store_everything, num_mcsim),
    m_baselines_orig(baselines) {

  m_total_count = arma::accu(counts);
  m_counts = arma::cumsum(counts);
  m_baselines = arma::cumsum(baselines);

  store = (store_everything ? &PBPOIscan::store_all : &PBPOIscan::store_max);

  // Reserve sizes for values calculated on observed data
  m_relrisk_in.set_size(m_out_length);
  m_relrisk_out.set_size(m_out_length);

  // Reserve sizes for values calculated on simulated
  sim_relrisk_in.set_size(m_num_mcsim);
  sim_relrisk_out.set_size(m_num_mcsim);
}

// Workhorse functions ---------------------------------------------------------

inline void PBPOIscan::calculate(const int storage_index,
                                 const int zone_nr,
                                 const int duration,
                                 const arma::uvec& current_zone,
                                 const arma::uvec& current_rows) {
  double C, B, risk_in, risk_out, term2, score;

  arma::uvec row_idx = current_rows.tail(1);

  // Counts and baselines are already aggregated
  C = arma::accu(m_counts.submat(row_idx, current_zone));
  B = arma::accu(m_baselines.submat(row_idx, current_zone));

  risk_in  = C / B;
  risk_out = (m_total_count > B ?
             (m_total_count - C) / (m_total_count - B) :
             1.0);
  term2 = (std::abs(m_total_count - C) < 1e-16 ?
          0.0 : (m_total_count - C) * std::log(risk_out));

  score = C > B ? C * std::log(risk_in) + term2 : 0.0;

  (this->*store)(storage_index, score, risk_in, risk_out, zone_nr + 1,
                 duration + 1);
}

inline void PBPOIscan::simulate_counts() {
  Rcpp::NumericVector probs(m_counts.n_cols * m_counts.n_rows);
  probs = armaVec2rcppVec(arma::vectorise(m_baselines_orig)) / m_total_count;
  // probs = Rcpp::as<Rcpp::NumericVector>(
  //           Rcpp::wrap(arma::vectorise(m_baselines_orig))) / m_total_count;

  arma::uvec vec_counts(m_counts.n_cols * m_counts.n_rows);
  vec_counts = rcppIVec2armaVec(
    Rcpp::RcppArmadillo::Rf_rmultinom(m_total_count, probs));
  // vec_counts = Rcpp::as<arma::uvec>(
  //   Rcpp::RcppArmadillo::Rf_rmultinom(m_total_count, probs));

  for (arma::uword j = 0; j < m_counts.n_cols; ++j) {
    m_counts.col(j) = arma::cumsum(
      vec_counts.subvec(j * m_counts.n_rows, (j + 1) * m_counts.n_rows - 1));
  }
}

inline int PBPOIscan::draw_sample(arma::uword row, arma::uword col) {
  return 1;
}

// Storage functions -----------------------------------------------------------

inline void PBPOIscan::store_all(int storage_index, double score, double q_in,
                                 double q_out, int zone_nr, int duration) {
  m_scores[storage_index]       = score;
  m_relrisk_in[storage_index]   = q_in;
  m_relrisk_out[storage_index]  = q_out;
  m_zone_numbers[storage_index] = zone_nr;
  m_durations[storage_index]    = duration;
}

inline void PBPOIscan::store_max(int storage_index, double score, double q_in,
                                 double q_out, int zone_nr, int duration) {
  if (score > m_scores[0]) {
    m_scores[0]       = score;
    m_relrisk_in[0]   = q_in;
    m_relrisk_out[0]  = q_out;
    m_zone_numbers[0] = zone_nr;
    m_durations[0]    = duration;
  }
}

inline void PBPOIscan::store_sim(int storage_index, double score, double q_in,
                                 double q_out, int zone_nr, int duration) {
  if (score > sim_scores[m_mcsim_index]) {
    sim_scores[m_mcsim_index]       = score;
    sim_relrisk_in[m_mcsim_index]   = q_in;
    sim_relrisk_out[m_mcsim_index]   = q_out;
    sim_zone_numbers[m_mcsim_index] = zone_nr;
    sim_durations[m_mcsim_index]    = duration;
  }
}

inline void PBPOIscan::set_sim_store_fun() {
  store = &PBPOIscan::store_sim;
}

// Retrieval functions ---------------------------------------------------------

inline Rcpp::DataFrame PBPOIscan::get_scan() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")         = m_zone_numbers,
    Rcpp::Named("duration")     = m_durations,
    Rcpp::Named("score")        = m_scores,
    Rcpp::Named("relrisk_in")   = m_relrisk_in,
    Rcpp::Named("relrisk_out")  = m_relrisk_out);
}

inline Rcpp::DataFrame PBPOIscan::get_mcsim() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")         = sim_zone_numbers,
    Rcpp::Named("duration")     = sim_durations,
    Rcpp::Named("score")        = sim_scores,
    Rcpp::Named("relrisk_in")   = sim_relrisk_in,
    Rcpp::Named("relrisk_out")  = sim_relrisk_out);
}

#endif
