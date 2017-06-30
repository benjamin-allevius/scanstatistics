#ifndef PBPOISCAN_H
#define PBPOISCAN_H

#include "USTscan.h"

class PBPOIscan : public USTscan<arma::umat> {

public:
  PBPOIscan(const arma::umat& counts, // counts = apply(obs_counts, 2, cumsum)
            const arma::mat& baselines, // basel = apply(obs_basel, 2, cumsum)
            const int total_count,
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
  int m_total_count;

  // Components of returned list
  arma::uvec m_zone_numbers;
  arma::uvec m_durations;
  arma::vec  m_scores;
  arma::vec  m_relrisk_in;
  arma::vec  m_relrisk_out;

  // Functions
  using store_ptr = void (PBPOIscan::*)(int storage_index, double score,
                                        double q_in, double q_out, int zone_nr,
                                        int duration);
  store_ptr store;
  void store_max(int storage_index, double score, double q_in, double q_out,
                 int zone_nr, int duration);
  void store_all(int storage_index, double score, double q_in, double q_out,
                 int zone_nr, int duration);

};

// Implementations -------------------------------------------------------------

inline PBPOIscan::PBPOIscan(const arma::umat& counts,
                            const arma::mat& baselines,
                            const int total_count,
                            const arma::uvec& zones,
                            const arma::uvec& zone_lengths,
                            const int num_locs,
                            const int num_zones,
                            const int max_dur,
                            const bool store_everything)
  : USTscan(counts, zones, zone_lengths, num_locs, num_zones, max_dur),
    m_baselines(baselines),
    m_total_count(total_count) {

  int out_length;
  if (store_everything) {
    out_length = m_num_zones * m_max_dur;
    store = &PBPOIscan::store_all;
  } else {
    out_length = 1;
    store = &PBPOIscan::store_max;
  }

  m_zone_numbers.set_size(out_length);
  m_durations.set_size(out_length);
  m_scores.set_size(out_length);
  m_relrisk_in.set_size(out_length);
  m_relrisk_out.set_size(out_length);

  if (!store_everything) {
    m_scores[0] = -1.0;
  }
}

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

inline Rcpp::DataFrame PBPOIscan::get_results() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")         = m_zone_numbers,
    Rcpp::Named("duration")     = m_durations,
    Rcpp::Named("score")        = m_scores,
    Rcpp::Named("relrisk_in")   = m_relrisk_in,
    Rcpp::Named("relrisk_out")  = m_relrisk_out);
}

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

#endif
