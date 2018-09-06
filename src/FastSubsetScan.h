#ifndef FASTSUBSETSCAN_H
#define FASTSUBSETSCAN_H

#include <cmath>
#include <vector>
#include <utility>
#include "probability_functions.h"
#include "RcppArmadillo.h"
// [[depends(RcppArmadillo)]]

using locationSubset = arma::uvec;
using streamSubset = arma::uvec;

struct Subset {
  arma::uvec locations;
  arma::uvec streams;
  arma::uword duration;

  Subset(arma::uvec locs, arma::uvec strms, arma::uword dur) :
    locations {locs}, streams {strms}, duration {dur} {}
};

// Fast Subset Scan Base class -------------------------------------------------

template <T>
class FastSubsetScan {

public:
  FastSubsetScan(
    const arma::cube<T>& counts,
    const arma::uvec& zones,
    const arma::uvec& zone_lengths,
    const arma::uvec& streams,
    const arma::uvec& stream_lengths,
    const bool store_everything,
    const arma::uword num_mcsim);
  
  virtual void run_scan() = 0;
  virtual void run_mcsim() = 0;

  virtual Rcpp::List get_scan() = 0;
  virtual Rcpp::List get_mcsim() = 0;

protected:
  arma::cube<T> m_counts;
  
  arma::uvec    m_zones;
  arma::uvec    m_zone_lengths;
  arma::uvec    m_streams;
  arma::uvec    m_stream_lengths;
  
  bool          m_store_everything;
  arma::uword   m_num_mcsim;
  
  arma::uword   m_num_locs;
  arma::uword   m_num_zones;
  arma::uword   m_num_streams;
  arma::uword   m_max_dur;
  arma::uword   m_mcsim_index;

  // Values calculated on observed data
  std::vector<Subset> m_subsets;
  arma::vec  m_scores;

  // Values calculated on simulated data
  std::vector<Subset> sim_subsets;
  arma::vec  sim_scores;

  // Functions
  virtual double score_fun(Subset s) = 0;
  virtual double priority_fun(Subset s) = 0;
  
  // Algorithm 1: optimize over all location subsets (zones) for fixed stream 
  //              subset and duration
  virtual 
  std::pair<Subset, double> 
  optimize_locations_quickly(const streamSubset& streams, 
                             const arma::uword duration) = 0;
  
  // Algorithm 2: optimize over all stream subsets for fixed location subset and
  //              duration
  virtual
  std::pair<Subset, double> 
  optimize_streams_quickly(const locationSubset& locations, 
                           const arma::uword duration) = 0;
  
  // Algorithm 3: optimize over stream subsets naively (slow), but fast over
  //              location subsets
  virtual std::pair<Subset, double> optimize_streams_naively() = 0;
  
  // Algorithm 4: optimize over locations subsets naively (slow), but fast over
  //              stream subsets
  virtual std::pair<Subset, double> optimize_locations_naively() = 0;
  
  // Algorithm 5: optimize over both location and stream subsets quickly
  virtual std::pair<Subset, double> optimize_both_quickly() = 0;

};

template <class T>
inline FastSubsetScan<T>::FastSubsetScan(
  const arma::cube<T>& counts,
  const arma::uvec& zones,
  const arma::uvec& zone_lengths,
  const arma::uvec& streams,
  const arma::uvec& stream_lengths,
  const bool store_everything,
  const arma::uword num_mcsim) :
    m_counts(counts),
    m_zones(zones),
    m_zone_lengths(zone_lengths),
    m_streams(streams),
    m_stream_lengths(stream_lengths),
    m_store_everything(store_everything),
    m_num_mcsim(num_mcsim) 
    m_num_locs(counts.n_cols),
    m_num_zones(zone_lengths.n_elem),
    m_num_streams(counts.n_slices),
    m_max_dur(counts.n_rows),
    m_mcsim_index(0) {
}

// Subset Aggregation ----------------------------------------------------------

template <class T>
class SubsetAggregation : public FastSubsetScan<T> {
public:
  SubsetAggregation(
    const arma::cube<T>& counts,
    const arma::uvec& zones,
    const arma::uvec& zone_lengths,
    const arma::uvec& streams,
    const arma::uvec& stream_lengths,
    const bool store_everything,
    const arma::uword num_mcsim);

protected:
};

template <class T>
inline SubsetAggregation<T>::SubsetAggregation(
  const arma::cube<T>& counts,
  const arma::uvec& zones,
  const arma::uvec& zone_lengths,
  const arma::uvec& streams,
  const arma::uvec& stream_lengths,
  const bool store_everything,
  const arma::uword num_mcsim) : 
    FastSubsetScan<T>(counts, 
                      zones, 
                      zone_lengths,
                      streams,
                      stream_lengths
                      store_everything, 
                      num_mcsim) {

}

template <T>
inline std::pair<Subset, double> 
SubsetAggregation<T>::optimize_locations_quickly(const streamSubset& streams, 
                                                 const arma::uword duration) {
  // 1. Compute time and stream aggregates for each location.
  // 2. Compute location priorities.
  // 3. Sort locations according to priorities.
  // 4. For increasing subsets of priority-ordered locations, compute scores.
  // 5. Return top scoring subset and its score.
}

template <T>
inline std::pair<Subset, double> 
SubsetAggregation<T>::optimize_streams_quickly(const locationSubset& locations, 
                                               const arma::uword duration) {
  // 1. Compute time and location aggregates for each stream.
  // 2. Compute stream priorities.
  // 3. Sort streams according to priorities.
  // 4. For increasing subsets of priority-ordered streams, compute scores.
  // 5. Return top scoring subset and its score.
}

template <T>
inline std::pair<Subset, double> 
SubsetAggregation<T>::optimize_streams_naively();
  // For each stream subset and duration:
  //   1. Run Algorithm 1 (optimize_locations_quickly).
  //   2. Return top scoring subset and its score.
}

template <T>
inline std::pair<Subset, double> 
SubsetAggregation<T>::optimize_streams_naively();
// For each location subset and duration:
//   1. Run Algorithm 2 (optimize_streams_quickly).
//   2. Return top scoring subset and its score.
}

template <T>
inline std::pair<Subset, double> 
SubsetAggregation<T>::optimize_both_quickly();
// For each duration:
//   1. Randomly choose a stream subset.
//   2. Get top-scoring location subset S* using Algorithm 1.
//   3. Use Algorithm 2 with S*.
//   4. Alternate between steps 2 and 3 until the score converges.
//   5. Repeat steps 1-5 a number of times.
// Return the top-scoring subset over all durations.
}


// Score Aggregation -----------------------------------------------------------

template <class T>
class ScoreAggregation : public FastSubsetScan<T> {
public:
  ScoreAggregation(const arma::cube<T>& counts,
                   const arma::uvec& zones,
                   const arma::uvec& zone_lengths,
                   const bool store_everything,
                   const arma::uword num_mcsim);
  
protected:
  // Algorithm 1: optimize over all stream subsets for fixed stream subset and 
  //              duration
  std::pair<Subset, double> 
  optimize_streams_quickly(const locationSubset& streams, 
                           const arma::uword duration);
  
  // Algorithm 2: optimize over all stream subsets for fixed location subset and
  //              duration
  optimize_locations_quickly(const streamSubset& streams, 
                             const arma::uword duration);
  
  // Algorithm 3: optimize over all stream subsets for fixed location subset and
  //              duration
  std::pair<Subset, double> optimize_locations_naively();
  
  // Algorithm 4: optimize over both location and stream subsets quickly
  std::pair<Subset, double> optimize_both_quickly();
};

template <T>
inline std::pair<Subset, double> 
ScoreAggregation<T>::optimize_streams_quickly(const locationSubset& streams, 
                                              const arma::uword duration) {
  // 1. Compute location-time aggregates for each stream.
  // 2. Compute scores using these aggregates.
  // 3. Return the subset of streams that has positive scores.
}

template <T>
inline std::pair<Subset, double> 
ScoreAggregation<T>::optimize_locations_naively() {
  // 1. For each location subset (zone) and duration, run Algorithm 1.
  // 2. Return the top-scoring subset.
}

template <T>
inline std::pair<Subset, double> 
ScoreAggregation<T>::optimize_both_naively() {
  // For each duration W:
  //   1. Choose p ~ Unif(0, 1) and include each stream w.p. p (but include at 
  //      least one stream).
  //   2. Initialize q_m for included streams, e.g. by q_m = exp(x_m), 
  //      x_m ~ Unif(0, 2). Set q_m = 1 for non-included streams.
  //   3. Get optimal location subset S* by Algorithm 2.
  //   4. Estimate all q_m given S* and W.
  //   5. Repeat steps 3 and 4 until scores converge.
  //   6. Repeat steps 1-5 a number of times and return the top-scoring subset.
  // Return the top-scoring subset over all durations.
}


#endif
