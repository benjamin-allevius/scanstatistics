#include "scan_eb_poisson.h"
#include "EBPOIscan.h"

Rcpp::DataFrame scan_eb_poisson_cpp(const arma::umat& counts,
                                const arma::mat& baselines,
                                const arma::uvec& zones,
                                const arma::uvec& zone_lengths,
                                const int num_locs,
                                const int num_zones,
                                const int max_dur,
                                const bool store_everything) {

  EBPOIscan ob {counts, baselines, zones, zone_lengths, num_locs, num_zones,
                max_dur, store_everything};
  ob.run_scan();
  return ob.get_results();
}

