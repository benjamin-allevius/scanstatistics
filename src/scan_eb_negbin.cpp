#include "scan_eb_negbin.h"
#include "EBNBscan.h"

Rcpp::DataFrame scan_eb_negbin_cpp(const arma::umat& counts,
                                   const arma::mat& baselines,
                                   const arma::mat& overdisp,
                                   const arma::uvec& zones,
                                   const arma::uvec& zone_lengths,
                                   const int num_locs,
                                   const int num_zones,
                                   const int max_dur,
                                   const bool store_everything,
                                   const bool score_type) {

  EBNBscan ob {counts, baselines, overdisp, zones, zone_lengths, num_locs,
                num_zones, max_dur, store_everything, score_type};
  ob.run_scan();
  return ob.get_results();
}

