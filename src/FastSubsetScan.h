// #ifndef FASTSUBSETSCAN_H
// #define FASTSUBSETSCAN_H
// 
// #include <cmath>
// #include <vector>
// #include "probability_functions.h"
// #include "RcppArmadillo.h"
// // [[depends(RcppArmadillo)]]
// 
// template <class T, class t>
// class FastSubsetScan {
//   
// public:
//   FastSubsetScan(const T& counts,
//                  const arma::uvec& m_zones,
//                  const arma::uvec& m_zone_lengths,
//                  const bool store_everything,
//                  const int num_mcsim);
//   void run_scan() = 0;
//   virtual void run_mcsim() = 0;
//   
//   virtual Rcpp::List get_scan() = 0;
//   virtual Rcpp::List get_mcsim() = 0;
//   
// protected:
//   int        m_num_locs;
//   int        m_num_zones;
//   int        m_num_streams;
//   int        m_max_dur;
//   int        m_num_mcsim;
//   bool       m_store_everything;
//   int        m_mcsim_index;
//   int        m_out_length;
//   T          m_counts;
//   arma::cube m_baselines;
//   arma::uvec m_zones;
//   arma::uvec m_zone_lengths;
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
//   struct Subset {
//     arma::uvec locations;
//     arma::uvec streams;
//     int duration;
//     double score;
//     
//     Subset(arma::uvec loc, arma::uvec strms, int dur, double scr) :
//       locations {loc}, streams {strms}, duration {dur}, score {scr} {}
//   };
//   
// };
// 
// template <class T, class t>
// class SubsetAggregation : public FastSubsetScan<T, t> {
// public:
//   
// protected:
//   arma::uvec optimize_locations(const arma::uvec& streams, 
//                                 const arma::uword w);
//   arma::uvec optimize_streams(const arma::uvec& streams, 
//                               const arma::uword w);
//   Subset fastloc_naivestream();
//   Subset naiveloc_faststream();
//   std::vector<Subset> optimize_both();
// };
// 
// #endif
