#ifndef EBPOISCAN_H
#define EBPOISCAN_H

#include "USTscan.h"

class EBPOIscan : public USTscan<arma::umat, arma::uword> {

public:
  EBPOIscan(const arma::umat& counts,
            const arma::mat& baselines,
            const arma::uvec& zones,
            const arma::uvec& zone_lengths,
            const bool store_everything,
            const arma::uword num_mcsim);

  Rcpp::DataFrame get_scan()  override;
  Rcpp::DataFrame get_mcsim() override;

private:
  arma::mat m_baselines;
  arma::mat m_baselines_orig; // used for simulation

  // Values calculated on observed data
  arma::vec m_relrisks;

  // Values calculated on simulated data
  arma::vec sim_relrisks;

  // Functions
  void calculate(const arma::uword storage_index,
                 const arma::uword zone_nr,
                 const arma::uword duration,
                 const arma::uvec& current_zone,
                 const arma::uvec& current_rows) override;
  void simulate_counts() override;
  arma::uword draw_sample(arma::uword row, arma::uword col) override;
  void set_sim_store_fun() override;
  
  using store_ptr = void (EBPOIscan::*)(arma::uword storage_index, double score,
                                        double q, arma::uword zone_nr, 
                                        arma::uword duration);
  store_ptr store;
  void store_max(arma::uword storage_index, double score, double q, 
                 arma::uword zone_nr, arma::uword duration);
  void store_all(arma::uword storage_index, double score, double q,
                 arma::uword zone_nr, arma::uword duration);
  void store_sim(arma::uword storage_index, double score, double q, 
                 arma::uword zone_nr, arma::uword duration);

};

// Implementations -------------------------------------------------------------

inline EBPOIscan::EBPOIscan(const arma::umat& counts,
                            const arma::mat& baselines,
                            const arma::uvec& zones,
                            const arma::uvec& zone_lengths,
                            const bool store_everything,
                            const arma::uword num_mcsim)
  : USTscan(counts, zones, zone_lengths, store_everything, num_mcsim),
    m_baselines_orig(baselines) {
  
  m_counts = arma::cumsum(counts);
  m_baselines = arma::cumsum(baselines);

  store = (store_everything ? &EBPOIscan::store_all : &EBPOIscan::store_max);

  // Reserve sizes for values calculated on observed data
  m_relrisks.set_size(m_out_length);

  // Reserve sizes for values calculated on simulated
  sim_relrisks.set_size(m_num_mcsim);

}

// Workhorse functions ---------------------------------------------------------

inline void EBPOIscan::calculate(const arma::uword storage_index,
                                 const arma::uword zone_nr,
                                 const arma::uword duration,
                                 const arma::uvec& current_zone,
                                 const arma::uvec& current_rows) {

  arma::uvec row_idx = current_rows.tail(1);

  // Counts and baselines are already aggregated
  double C = arma::accu(m_counts.submat(row_idx, current_zone));
  double B = arma::accu(m_baselines.submat(row_idx, current_zone));

  double score = C > B ? C * log(C / B) + B - C : 0.0;

  (this->*store)(storage_index, score, std::max(1.0, C / B), zone_nr + 1,
   duration + 1);
}

inline arma::uword EBPOIscan::draw_sample(arma::uword row, arma::uword col) {
  return R::rpois(m_baselines_orig.at(row, col));
}

inline void EBPOIscan::simulate_counts() {
  for (arma::uword j = 0; j < m_counts.n_cols; ++j) {
    for (arma::uword i = 0; i < m_counts.n_rows; ++i) {
      m_counts.at(i, j) = draw_sample(i, j);
    }
  }
  m_counts = arma::cumsum(m_counts);
}

// Storage functions -----------------------------------------------------------

inline void EBPOIscan::store_all(arma::uword storage_index, double score, 
                                 double q, arma::uword zone_nr, 
                                 arma::uword duration) {
  m_scores[storage_index]       = score;
  m_relrisks[storage_index]     = q;
  m_zone_numbers[storage_index] = zone_nr;
  m_durations[storage_index]    = duration;
}

inline void EBPOIscan::store_max(arma::uword storage_index, double score, 
                                 double q, arma::uword zone_nr, 
                                 arma::uword duration) {
  if (score > m_scores[0]) {
    m_scores[0]       = score;
    m_relrisks[0]     = q;
    m_zone_numbers[0] = zone_nr;
    m_durations[0]    = duration;
  }
}

inline void EBPOIscan::store_sim(arma::uword storage_index, double score, 
                                 double q, arma::uword zone_nr, 
                                 arma::uword duration) {
  if (score > sim_scores[m_mcsim_index]) {
    sim_scores[m_mcsim_index]       = score;
    sim_relrisks[m_mcsim_index]     = q;
    sim_zone_numbers[m_mcsim_index] = zone_nr;
    sim_durations[m_mcsim_index]    = duration;
  }
}

inline void EBPOIscan::set_sim_store_fun() {
  store = &EBPOIscan::store_sim;
}

// Retrieval functions ---------------------------------------------------------

inline Rcpp::DataFrame EBPOIscan::get_scan() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")     = m_zone_numbers,
    Rcpp::Named("duration") = m_durations,
    Rcpp::Named("score")    = m_scores,
    Rcpp::Named("relrisk")  = m_relrisks);
}

inline Rcpp::DataFrame EBPOIscan::get_mcsim() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")     = sim_zone_numbers,
    Rcpp::Named("duration") = sim_durations,
    Rcpp::Named("score")    = sim_scores,
    Rcpp::Named("relrisk")  = sim_relrisks);
}



#endif
