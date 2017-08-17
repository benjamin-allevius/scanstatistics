// #ifndef BGPSCAN_H
// #define BGPSCAN_H
// 
// #include "USTscan.h"
// 
// class BGPscan : public USTscanbase<arma::umat, arma::uword> {
// 
// public:
//   BGPscan(const arma::umat& counts,
//           const arma::mat& baselines,
//           const arma::uvec& zones,
//           const arma::uvec& zone_lengths,
//           const bool store_everything);
// 
//   Rcpp::DataFrame get_scan() override;
// 
// private:
//   arma::mat m_baselines;
//   arma::mat m_baselines_orig; // used for simulation
//   arma::uword m_total_count;
// 
//   // Functions
//   void calculate(const arma::uword storage_index,
//                  const arma::uword zone_nr,
//                  const arma::uword duration,
//                  const arma::uvec& current_zone,
//                  const arma::uvec& current_rows) override;
// 
//   using store_ptr = void (BGPscan::*)(arma::uword storage_index, double score,
//                           arma::uword zone_nr, arma::uword duration);
//   store_ptr store;
//   void store_all(arma::uword storage_index, double score, arma::uword zone_nr,
//                  arma::uword duration);
//   
//   void post_process() override;
// 
// };
// 
// // Implementations -------------------------------------------------------------
// 
// inline BGPscan::BGPscan(const arma::umat& counts,
//                         const arma::mat& baselines,
//                         const arma::uvec& zones,
//                         const arma::uvec& zone_lengths,
//                         const bool store_everything)
//   : USTscanbase(counts, zones, zone_lengths, true),
//     m_baselines_orig(baselines) {
// 
//   m_total_count = arma::accu(counts);
//   m_counts = arma::cumsum(counts);
//   m_baselines = arma::cumsum(baselines);
// 
//   store = &BGPscan::store_all;
// 
// }
// 
// // Workhorse functions ---------------------------------------------------------
// 
// inline void BGPscan::calculate(const arma::uword storage_index,
//                                const arma::uword zone_nr,
//                                const arma::uword duration,
//                                const arma::uvec& current_zone,
//                                const arma::uvec& current_rows) {
// 
//   arma::uword C;
//   double B, score;
//   
//   arma::uvec row_idx = current_rows.tail(1);
//   
//   // Counts and baselines are already aggregated
//   C = arma::accu(m_counts.submat(row_idx, current_zone));
//   B = arma::accu(m_baselines.submat(row_idx, current_zone));
//   
//   
//   
//   (this->*store)(storage_index, score, zone_nr + 1, duration + 1);
// }
// 
// inline void BGPscan::post_process() {
//   
// }
// 
// 
// // Storage functions -----------------------------------------------------------
// 
// inline void BGPscan::store_all(arma::uword storage_index, double score,
//                                  arma::uword zone_nr, arma::uword duration) {
//   m_scores[storage_index]       = score;
//   m_zone_numbers[storage_index] = zone_nr;
//   m_durations[storage_index]    = duration;
// }
// 
// // Retrieval functions ---------------------------------------------------------
// 
// inline Rcpp::DataFrame BGPscan::get_scan() {
//   return Rcpp::DataFrame::create(
//     Rcpp::Named("zone")     = m_zone_numbers,
//     Rcpp::Named("duration") = m_durations,
//     Rcpp::Named("score")    = m_scores);
// }
// 
// #endif
