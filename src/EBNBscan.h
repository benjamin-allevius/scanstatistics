#ifndef EBNBSCAN_H
#define EBNBSCAN_H

#include "USTscan.h"
#include "scan_utility.h"

class EBNBscan : public USTscan<arma::umat, arma::uword> {

public:
  EBNBscan(const arma::umat& counts,
            const arma::mat& baselines,
            const arma::mat& overdisp,
            const arma::uvec& zones,
            const arma::uvec& zone_lengths,
            const bool store_everything,
            const arma::uword num_mcsim,
            const bool score_hotspot);

  Rcpp::DataFrame get_scan()  override;
  Rcpp::DataFrame get_mcsim() override;

private:
  arma::mat m_baselines;
  arma::mat m_overdisp;

  // Functions
  void calculate(const arma::uword storage_index,
                 const arma::uword zone_nr,
                 const arma::uword duration,
                 const arma::uvec& current_zone,
                 const arma::uvec& current_rows) override;

  arma::uword draw_sample(arma::uword row, arma::uword col) override;
  void set_sim_store_fun() override;

  using store_ptr = void (EBNBscan::*)(arma::uword storage_index, double score,
                                       arma::uword zone_nr, arma::uword duration);
  store_ptr store;
  void store_max(arma::uword storage_index, double score, arma::uword zone_nr, arma::uword duration);
  void store_all(arma::uword storage_index, double score, arma::uword zone_nr, arma::uword duration);
  void store_sim(arma::uword storage_index, double score, arma::uword zone_nr, arma::uword duration);

  using score_ptr = double (EBNBscan::*)(const arma::uvec& y,
                                         const arma::vec& mu,
                                         const arma::vec& omega,
                                         const arma::uword d);
  score_ptr score_fun;
  double score_emerge(const arma::uvec& y, const arma::vec& mu,
                      const arma::vec& omega, const arma::uword d);
  double score_hotspot(const arma::uvec& y, const arma::vec& mu,
                       const arma::vec& omega, const arma::uword d);

};

// Implementations -------------------------------------------------------------

inline EBNBscan::EBNBscan(const arma::umat& counts,
                          const arma::mat& baselines,
                          const arma::mat& overdisp,
                          const arma::uvec& zones,
                          const arma::uvec& zone_lengths,
                          const bool store_everything,
                          const arma::uword num_mcsim,
                          const bool score_hotspot)
  : USTscan(counts, zones, zone_lengths, store_everything, num_mcsim),
    m_baselines(baselines),
    m_overdisp(overdisp) {

  store = (m_store_everything ? &EBNBscan::store_all : &EBNBscan::store_max);

  score_fun = (score_hotspot ?
               &EBNBscan::score_hotspot : &EBNBscan::score_emerge);
}

// Workhorse functions ---------------------------------------------------------

inline void EBNBscan::calculate(const arma::uword storage_index,
                                const arma::uword zone_nr,
                                const arma::uword duration,
                                const arma::uvec& current_zone,
                                const arma::uvec& current_rows) {
  // Extract counts and parameters as vectors
  arma::uvec y = arma::vectorise(m_counts.submat(current_rows,
                                                 current_zone));
  arma::vec mu = arma::vectorise(m_baselines.submat(current_rows,
                                                    current_zone));
  arma::vec omega = arma::vectorise(m_overdisp.submat(current_rows,
                                                      current_zone));
  double score = (this->*score_fun)(y, mu, omega, duration);

  (this->*store)(storage_index, score, zone_nr + 1, duration + 1);
}

inline double EBNBscan::score_hotspot(const arma::uvec& y, const arma::vec& mu,
                                      const arma::vec& omega, const arma::uword d) {
  return arma::accu((y - mu) / omega) / arma::accu(mu / omega);
}

inline double EBNBscan::score_emerge(const arma::uvec& y, const arma::vec& mu,
                                     const arma::vec& omega, const arma::uword d) {
  double num = 0.0;
  double den = 0.0;
  arma::uword idx = 0;
  for (arma::uword i = 0; i < y.n_elem / (d + 1); ++i) {
    for (arma::uword t = 0; t < d + 1; ++idx, ++t) {
      num += (y[idx] - mu[idx]) * (d + 1.0 - t) / omega[idx];
      den += mu[idx] * std::pow(d + 1.0 - t, 2) / omega[idx];
    }
  }
  return num / den;
}

inline arma::uword EBNBscan::draw_sample(arma::uword row, arma::uword col) {
  return rnbinom2(m_baselines.at(row, col), m_overdisp.at(row, col));
}


// Storage functions -----------------------------------------------------------

inline void EBNBscan::store_all(arma::uword storage_index, double score, arma::uword zone_nr,
                                arma::uword duration) {
  m_scores[storage_index]       = score;
  m_zone_numbers[storage_index] = zone_nr;
  m_durations[storage_index]    = duration;
}

inline void EBNBscan::store_max(arma::uword storage_index, double score, arma::uword zone_nr,
                                arma::uword duration) {
  if (score > m_scores[0]) {
    m_scores[0]       = score;
    m_zone_numbers[0] = zone_nr;
    m_durations[0]    = duration;
  }
}

inline void EBNBscan::store_sim(arma::uword storage_index, double score, arma::uword zone_nr,
                                arma::uword duration) {
  if (score > sim_scores[m_mcsim_index]) {
    sim_scores[m_mcsim_index]       = score;
    sim_zone_numbers[m_mcsim_index] = zone_nr;
    sim_durations[m_mcsim_index]    = duration;
  }
}

inline void EBNBscan::set_sim_store_fun() {
  store = &EBNBscan::store_sim;
}

// Retrieval functions ---------------------------------------------------------

inline Rcpp::DataFrame EBNBscan::get_scan() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")     = m_zone_numbers,
    Rcpp::Named("duration") = m_durations,
    Rcpp::Named("score")    = m_scores);
}

inline Rcpp::DataFrame EBNBscan::get_mcsim() {
  return Rcpp::DataFrame::create(
    Rcpp::Named("zone")     = sim_zone_numbers,
    Rcpp::Named("duration") = sim_durations,
    Rcpp::Named("score")    = sim_scores);
}


#endif
