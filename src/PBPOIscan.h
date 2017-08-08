#ifndef PBPOISCAN_H
#define PBPOISCAN_H

#include "PBPOIabstract.h"
#include "scan_utility.h"
#include <RcppArmadilloExtensions/rmultinom.h>

class PBPOIscan : public PBPOIabstract {

public:
  PBPOIscan(const arma::umat& counts,
            const arma::mat& baselines,
            const arma::uvec& zones,
            const arma::uvec& zone_lengths,
            const bool store_everything,
            const arma::uword num_mcsim);

private:
  void simulate_counts() override;

};

// Implementations -------------------------------------------------------------

inline PBPOIscan::PBPOIscan(const arma::umat& counts,
                            const arma::mat& baselines,
                            const arma::uvec& zones,
                            const arma::uvec& zone_lengths,
                            const bool store_everything,
                            const arma::uword num_mcsim)
  : PBPOIabstract(counts, baselines, zones, zone_lengths, store_everything, 
                  num_mcsim) {}

// Workhorse functions ---------------------------------------------------------

inline void PBPOIscan::simulate_counts() {
  Rcpp::NumericVector probs(m_counts.n_cols * m_counts.n_rows);
  probs = vec2NumericVector(arma::vectorise(m_baselines_orig)) / m_total_count;

  arma::uvec vec_counts(m_counts.n_cols * m_counts.n_rows);
  vec_counts = IntegerVector2uvec(
    Rcpp::RcppArmadillo::Rf_rmultinom(m_total_count, probs));

  // Columns of m_counts should be cumulative sums
  for (arma::uword j = 0; j < m_counts.n_cols; ++j) {
    m_counts.col(j) = arma::cumsum(
      vec_counts.subvec(j * m_counts.n_rows, (j + 1) * m_counts.n_rows - 1));
  }
}

#endif
