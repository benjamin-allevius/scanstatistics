// #ifndef FASTSUBSETSCAN_H
// #define FASTSUBSETSCAN_H
// 
// #include <cmath>
// #include <vector>
// #include "probability_functions.h"
// #include "RcppArmadillo.h"
// // [[depends(RcppArmadillo)]]
// 
// struct Subset {
//   arma::uvec locations;
//   arma::uvec streams;
//   arma::uword duration;
//   double score;
// 
//   Subset(arma::uvec locs, arma::uvec strs, arma::uword dur, double scr) :
//     locations {locs}, streams {strs}, duration {dur}, score {scr} {}
// };
// 
// // CT = Cube type, MT = Matrix type, VT = vector type, T = scalar type
// template <class CT, class MT, class VT, class T>
// class FastSubsetScan {
// 
// public:
//   FastSubsetScan(const CT& counts,
//                  const arma::cube baselines,
//                  const arma::uvec& zones,
//                  const arma::uvec& zone_lengths,
//                  const bool store_everything,
//                  const arma::uword num_mcsim);
//   void run_scan() = 0;
//   virtual void run_mcsim() = 0;
// 
//   virtual Rcpp::List get_scan() = 0;
//   virtual Rcpp::List get_mcsim() = 0;
// 
// protected:
//   arma::uword m_num_locs;
//   arma::uword m_num_zones;
//   arma::uword m_num_streams;
//   arma::uword m_max_dur;
//   arma::uword m_num_mcsim;
//   bool        m_store_everything;
//   arma::uword m_mcsim_index;
//   arma::uword m_out_length;
//   CT    m_counts;
//   arma::cube  m_baselines;
//   arma::uvec  m_zones;
//   arma::uvec  m_zone_lengths;
// 
//   // Values calculated on observed data
//   arma::uvec m_zone_numbers;
//   arma::uvec m_durations;
//   arma::uvec m_streams;
//   arma::vec  m_scores;
// 
//   // Values calculated on simulated data
//   arma::uvec sim_zone_numbers;
//   arma::uvec sim_durations;
//   arma::uvec sim_streams;
//   arma::vec  sim_scores;
// 
//   // Functions
//   virtual double score_fun(T c, double b) = 0;
// 
//   // struct Subset {
//   //   arma::uvec locations;
//   //   arma::uvec streams;
//   //   int duration;
//   //   double score;
//   //
//   //   Subset(arma::uvec loc, arma::uvec strms, int dur, double scr) :
//   //     locations {loc}, streams {strms}, duration {dur}, score {scr} {}
//   // };
// 
// };
// 
// template <class CT, class MT, class VT, class T>
// inline FastSubsetScan<CT, MT, VT, T>::FastSubsetScan(
//     const CT& counts,
//     const arma::cube baselines,
//     const arma::uvec& zones,
//     const arma::uvec& zone_lengths,
//     const bool store_everything,
//     const arma::uword num_mcsim) :
//       m_counts(counts),  
//       m_num_locs(counts.n_cols),
//       m_num_zones(zone_lengths.n_elem),
//       m_num_streams(counts.n_slices),
//       m_max_dur(counts.n_rows),
//       m_zones(zones),
//       m_zone_lengths(zone_lengths),
//       m_store_everything(store_everything),
//       m_num_mcsim(num_mcsim) {
//   
//   
//   
// }
// 
// 
// template <class CT, class MT, class VT, class T>
// class SubsetAggregation : public FastSubsetScan<CT, MT, VT, T> {
// public:
//   SubsetAggregation(const CT& counts,
//                     const arma::cube baselines,
//                     const arma::uvec& zones,
//                     const arma::uvec& zone_lengths,
//                     const bool store_everything,
//                     const arma::uword num_mcsim);
// 
// protected:
//   arma::uvec opt_locs(const arma::uvec& streams, const arma::uword w);
//   arma::uvec opt_strs(const arma::uvec& streams, const arma::uword w);
//   Subset fastloc_naivestream();
//   Subset naiveloc_faststream();
//   std::vector<Subset> optimize_both();
// };
// 
// template <class CT, class MT, class VT, class T>
// inline SubsetAggregation<CT, MT, VT, T>::SubsetAggregation(
//       const CT& counts,
//       const arma::cube baselines,
//       const arma::uvec& zones,
//       const arma::uvec& zone_lengths,
//       const bool store_everything,
//       const arma::uword num_mcsim)
//   : FastSubsetScan<CT, MT, VT, T>(counts, baselines, zones, zone_lengths, 
//                                          store_everything, num_mcsim) {
//   
// }
// 
// template <class CT, class MT, class VT, class T>
// inline arma::uvec SubsetAggregation<CT, MT, VT, T>::opt_locs(
//     const arma::uvec& streams, const arma::uword w) {
//   VT C_i(m_num;
//   VT B_i;
// }
// 
// 
// #endif
