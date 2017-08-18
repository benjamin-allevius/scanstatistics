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
//           const bool store_everything,
//           const double outbreak_prob,
//           const double alpha_null,
//           const double beta_null,
//           const double alpha_alt,
//           const double beta_alt);
// 
//   Rcpp::DataFrame get_scan() override;
// 
// private:
//   arma::mat m_baselines;
//   arma::mat m_baselines_orig;
//   arma::uword m_total_count;
//   double m_total_baseline;
//   
//   // Probabilities
//   double m_log_outbreak_prob;
//   double m_log_null_prob;
//   double m_log_window_prob;
//   double m_log_data_prob;
//   
//   // Relative risk gamma distribution parameters
//   double m_alpha_null;
//   double m_beta_null;
//   double m_alpha_alt;
//   double m_beta_alt;
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
//   double log_prob(const arma::uword C, const double B, const double alpha, 
//                   const double beta);
// 
// };
// 
// // Implementations -------------------------------------------------------------
// 
// inline BGPscan::BGPscan(const arma::umat& counts,
//                         const arma::mat& baselines,
//                         const arma::uvec& zones,
//                         const arma::uvec& zone_lengths,
//                         const bool store_everything,
//                         const double outbreak_prob,
//                         const double alpha_null,
//                         const double beta_null,
//                         const double alpha_alt,
//                         const double beta_alt)
//   : USTscanbase(counts, zones, zone_lengths, true),
//     m_baselines_orig(baselines),
//     m_log_outbreak_prob(std::log(outbreak_prob)),
//     m_log_null_prob(std::log(1.0 - outbreak_prob)),
//     m_alpha_null(alpha_null),
//     m_beta_null(beta_null),
//     m_alpha_alt(alpha_alt),
//     m_beta_alt(beta_alt) {
// 
//   m_total_count = arma::accu(counts);
//   m_total_baseline = arma::accu(baselines);
//   
//   m_counts = arma::cumsum(counts);
//   m_baselines = arma::cumsum(baselines);
//   
//   m_log_window_prob = m_log_outbreak_prob 
//                     - std::log(m_num_zones) 
//                     - std::log(m_max_dur);
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
//   score = log_prob(C, B, m_alpha_alt, m_beta_alt) + 
//           log_prob(m_total_count - C, m_total_baseline - B, 
//                    m_alpha_null, m_beta_null) + m_log_window_prob;
// 
//   (this->*store)(storage_index, score, zone_nr + 1, duration + 1);
// }
// 
// inline void BGPscan::post_process() {
//   double top_score = m_scores.max();
//   double null_log_prob = log_prob(m_total_count, m_total_baseline,
//                                   m_alpha_null, m_beta_null) 
//                        + m_log_null_prob;
//   
//   // Log-Sum-Exp
//   double overall_max = std::max(top_score, null_log_prob);
//   double log_data_prob = std::exp(null_log_prob - overall_max);
//   for (arma::uword i = 0; i < m_scores.n_elem; ++i) {
//     log_data_prob += std::exp(m_scores.at(i) - overall_max);
//   }
//   log_data_prob = overall_max + std::log(log_data_prob);
//   m_log_data_prob = log_data_prob;
// }
// 
// inline double BGPscan::log_prob(const arma::uword C, const double B, 
//                                 const double alpha, const double beta) {
//   return alpha * std::log(beta) * std::lgamma(alpha + C) -
//          (alpha + C) * std::log(beta + B) - std::lgamma(alpha);
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
