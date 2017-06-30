#include "scan_eb_zip.h"
#include "EBZIPscan.h"

Rcpp::DataFrame scan_eb_zip_cpp(const arma::umat& counts,
                                const arma::mat& baselines,
                                const arma::mat& probs,
                                const arma::uvec& zones,
                                const arma::uvec& zone_lengths,
                                const int num_locs,
                                const int num_zones,
                                const int max_dur,
                                const double rel_tol,
                                const bool store_everything) {

  EBZIPscan ob {counts, baselines, probs, zones, zone_lengths, num_locs,
                num_zones, max_dur, rel_tol, store_everything};
  ob.run_scan();
  return ob.get_results();
}

