#ifndef EBPOISCAN_H
#define EBPOISCAN_H

#include "USTscan.h"

class EBPOIscan : public USTscan<arma::umat> {

public:
  EBPOIscan(const arma::umat& counts,
            const arma::mat& baselines,
            const arma::uvec& zones,
            const arma::uvec& zone_lengths,
            const int num_locs,
            const int num_zones,
            const int max_dur,
            const bool store_everything);
  void calculate(const int storage_index,
                 const int zone_nr,
                 const int duration,
                 const arma::uvec& current_zone,
                 const arma::uvec& current_rows);
  Rcpp::DataFrame get_results();

private:
  arma::mat m_baselines;

  // Components of returned list
  arma::uvec m_zone_numbers;
  arma::uvec m_durations;
  arma::vec  m_scores;
  arma::vec  m_relrisks;

  // Functions
  using store_ptr = void (EBPOIscan::*)(int storage_index, double score,
                                        double q, int zone_nr, int duration);
  store_ptr store;
  void store_max(int storage_index, double score, double q, int zone_nr,
                 int duration);
  void store_all(int storage_index, double score, double q, int zone_nr,
                 int duration);

};

// Implementations -------------------------------------------------------------

inline EBPOIscan::EBPOIscan(const arma::umat& counts,
                            const arma::mat& baselines,
                            const arma::uvec& zones,
                            const arma::uvec& zone_lengths,
                            const int num_locs,
                            const int num_zones,
                            const int max_dur,
                            const bool store_everything)
  : USTscan(counts, zones, zone_lengths, num_locs, num_zones, max_dur),
    m_baselines(baselines) {

  int out_length;
  if (store_everything) {
    out_length = m_num_zones * m_max_dur;
    store = &EBPOIscan::store_all;
  } else {
    out_length = 1;
    store = &EBPOIscan::store_max;
  }

  m_zone_numbers.set_size(out_length);
  m_durations.set_size(out_length);
  m_scores.set_size(out_length);
  m_relrisks.set_size(out_length);

  if (!store_everything) {
    m_scores[0] = -1.0;
  }
}

inline void EBPOIscan::calculate(const int storage_index,
                                 const int zone_nr,
                                 const int duration,
                                 const arma::uvec& current_zone,
                                 const arma::uvec& current_rows) {

  arma::uvec row_idx = current_rows.tail(1);

  // Counts and baselines are already aggregated
  double C = arma::accu(m_counts.submat(row_idx, current_zone));
  double B = arma::accu(m_baselines.submat(row_idx, current_zone));

  double score = (C > B ? C * (log(C) - log(B)) + B - C : 0.0);

  (this->*store)(storage_index, score, std::max(1.0, C / B), zone_nr + 1,
   duration + 1);
}

inline Rcpp::DataFrame EBPOIscan::get_results() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")     = m_zone_numbers,
    Rcpp::Named("duration") = m_durations,
    Rcpp::Named("score")    = m_scores,
    Rcpp::Named("relrisk")  = m_relrisks);
}

inline void EBPOIscan::store_all(int storage_index, double score, double q,
                                 int zone_nr, int duration) {
  m_scores[storage_index]       = score;
  m_relrisks[storage_index]     = q;
  m_zone_numbers[storage_index] = zone_nr;
  m_durations[storage_index]    = duration;
}

inline void EBPOIscan::store_max(int storage_index, double score, double q,
                                 int zone_nr, int duration) {
  if (score > m_scores[0]) {
    m_scores[0]       = score;
    m_relrisks[0]     = q;
    m_zone_numbers[0] = zone_nr;
    m_durations[0]    = duration;
  }
}

#endif
