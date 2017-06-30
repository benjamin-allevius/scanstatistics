#include "scan_pb_poisson.h"
#include "PBPOIscan.h"

Rcpp::DataFrame scan_pb_poisson_cpp(const arma::umat& counts,
                                const arma::mat& baselines,
                                const int total_count,
                                const arma::uvec& zones,
                                const arma::uvec& zone_lengths,
                                const int num_locs,
                                const int num_zones,
                                const int max_dur,
                                const bool store_everything) {

  PBPOIscan ob {counts, baselines, total_count, zones, zone_lengths, num_locs,
                num_zones, max_dur, store_everything};
  ob.run_scan();
  return ob.get_results();
}

