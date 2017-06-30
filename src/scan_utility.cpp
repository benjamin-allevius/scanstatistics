#include "scan_utility.h"

std::vector<int> get_zero_indices(arma::uvec v) {
  std::vector<int> zero_idx;
  for (int i = 0; i < v.n_elem; ++i) {
    if (v[i] == 0) zero_idx.push_back(i);
  }
  return zero_idx;
}
