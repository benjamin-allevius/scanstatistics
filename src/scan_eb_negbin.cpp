#include "scan_eb_negbin.h"
#include "EBNBscan.h"

Rcpp::List scan_eb_negbin_cpp(const arma::umat& counts,
                              const arma::mat& baselines,
                              const arma::mat& overdisp,
                              const arma::uvec& zones,
                              const arma::uvec& zone_lengths,
                              const bool store_everything,
                              const arma::uword num_mcsim,
                              const bool score_hotspot) {

  EBNBscan ob {counts, baselines, overdisp, zones, zone_lengths, 
               store_everything, num_mcsim, score_hotspot};
  ob.run_scan();
  ob.run_mcsim();
  return Rcpp::List::create(
    Rcpp::Named("observed") = ob.get_scan(),
    Rcpp::Named("simulated") = ob.get_mcsim());
}

