#ifndef USTSCAN_H
#define USTSCAN_H

#include <cmath>
#include "probability_functions.h"
#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]

// Univariate Spate-Time scan statistic (base class) ---------------------------
template <class T>
class USTscan {

public:
  USTscan(const T& counts,
          const arma::uvec& m_zones,
          const arma::uvec& m_zone_lengths,
          const int num_locs,
          const int num_zones,
          const int max_dur);
  void run_scan();

protected:
  int        m_num_locs;
  int        m_num_zones;
  int        m_max_dur;
  T          m_counts;
  arma::uvec m_zones;
  arma::uvec m_zone_lengths;

  virtual void calculate(const int storage_index,
                         const int zone_nr,
                         const int duration,
                         const arma::uvec& current_zone,
                         const arma::uvec& current_rows) = 0;

  arma::uvec zone_complement(const arma::uvec& zone);
  arma::uvec duration_complement(const int d);

};

// Implementations -------------------------------------------------------------

template <class T>
inline USTscan<T>::USTscan(const T& counts,
                           const arma::uvec& zones,
                           const arma::uvec& zone_lengths,
                           const int num_locs,
                           const int num_zones,
                           const int max_dur)
  : m_counts(counts),
    m_num_locs(num_locs),
    m_num_zones(num_zones),
    m_max_dur(max_dur),
    m_zones(zones),
    m_zone_lengths(zone_lengths) {

  // Do something here?
}


template <class T>
inline void USTscan<T>::run_scan() {
  int i = 0; // Storage index
  for (int d = 0; d < m_max_dur; ++d) {

    arma::uvec row_idx (d + 1); // d+1 latest time periods
    for (int k = 0; k <= d; ++k) row_idx[k] = k;

    // Indices for extracting the current zone
    int zone_start = 0;
    int zone_end = 0;

    for (int z = 0; z < m_num_zones; ++z) {

      // Extract zone
      zone_end = zone_start + m_zone_lengths[z] - 1;
      arma::uvec current_zone = m_zones(arma::span(zone_start, zone_end));

      // Calculate scan statistic + related quantities for current ST-window
      calculate(i, z, d, current_zone, row_idx);

      if (i % 500 == 0) Rcpp::checkUserInterrupt();

      zone_start = zone_end + 1;
      ++i;
    }
  }
}

template <class T>
inline arma::uvec USTscan<T>::zone_complement(const arma::uvec& zone) {
  // If there is only a single zone
  if (zone.n_elem == m_num_locs) {
    return arma::regspace<arma::uvec>(1, 1, 0);
  }
  arma::uvec res(m_num_locs - zone.n_elem);
  int idx = 0;

  for (int i = 0; i < m_num_locs; ++i) {
    bool outside_zone = true;

    for (const int& j : zone) {
      if (i == j) {
        outside_zone = false;
        break;
      }
    }

    if (outside_zone) {
      res[idx] = i;
      ++idx;
    }
  }

  return res;
}


template <class T>
inline arma::uvec USTscan<T>::duration_complement(const int d) {
  if (d == m_max_dur - 1) return arma::regspace<arma::uvec>(1, 1, 0);

  arma::uvec row_idx(m_max_dur - (d + 1)); // d+1 latest time periods
  for (int i = 0, k = d + 1; k < m_max_dur; ++i, ++k) {
    row_idx[i] = k;
  }
  return row_idx;
}

#endif
