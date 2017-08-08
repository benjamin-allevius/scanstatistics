#ifndef PBPERMSCAN_H
#define PBPERMSCAN_H

#include "PBPOIabstract.h"
#include "scan_utility.h"

class PBPERMscan : public PBPOIabstract {

public:
  PBPERMscan(const arma::umat& counts,
             const arma::mat& baselines,
             const arma::uvec& zones,
             const arma::uvec& zone_lengths,
             const bool store_everything,
             const arma::uword num_mcsim);

private:
  // Each case in counts expanded --> use for permutation
  arma::umat m_counts_expanded;
  arma::uvec m_time_counts;
  arma::uvec m_loc_counts;
  void simulate_counts() override;

};

// Implementations -------------------------------------------------------------

inline PBPERMscan::PBPERMscan(const arma::umat& counts,
                              const arma::mat& baselines,
                              const arma::uvec& zones,
                              const arma::uvec& zone_lengths,
                              const bool store_everything,
                              const arma::uword num_mcsim)
  : PBPOIabstract(counts, baselines, zones, zone_lengths, store_everything, 
                  num_mcsim) {
  m_counts_expanded = expand_matrix(counts);
  m_time_counts = m_counts_expanded.col(0);
  // m_loc_counts = m_counts_expanded.col(1);
  
}

// Workhorse functions ---------------------------------------------------------

inline void PBPERMscan::simulate_counts() {
  m_counts_expanded.col(0) = shuffle_time_counts(m_time_counts);
  m_counts = arma::cumsum(contract_matrix(m_counts_expanded,
                                          m_counts.n_rows, m_counts.n_cols));
}

#endif
