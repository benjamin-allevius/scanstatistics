#include "scan_pb_perm.h"
#include "PBPERMscan.h"

Rcpp::List scan_pb_perm_cpp(const arma::umat& counts,
                            const arma::mat& baselines,
                            const arma::uvec& zones,
                            const arma::uvec& zone_lengths,
                            const bool store_everything,
                            const arma::uword num_mcsim) {

  PBPERMscan ob {counts, baselines, zones, zone_lengths, store_everything, 
                 num_mcsim};
  ob.run_scan();
  ob.run_mcsim();
  return Rcpp::List::create(
    Rcpp::Named("observed") = ob.get_scan(),
    Rcpp::Named("simulated") = ob.get_mcsim());
}

