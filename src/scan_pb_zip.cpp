#include "scan_pb_zip.h"
#include "PBZIPscan.h"

Rcpp::DataFrame scan_pb_zip_cpp(const arma::umat& counts,
                                const arma::mat& pop,
                                const arma::uvec& zones,
                                const arma::uvec& zone_lengths,
                                const int num_locs,
                                const int num_zones,
                                const int max_dur,
                                const double rel_tol,
                                const bool store_everything) {

  PBZIPscan ob {counts, pop, zones, zone_lengths, num_locs,
                num_zones, max_dur, rel_tol, store_everything};
  ob.run_scan();
  return ob.get_results();
}

