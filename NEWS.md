# scanstatistics 0.2

# Major changes

* Multiple functions reimplemented in C++.
* Fast detection for multivariate space-time data with the function 
  `subset_aggregation`, for multiple distributional assumptions.

## Minor changes

* The functions `knn_zones` and `flexible_zones` now run faster due to change
  in algorithms.
* Bug fix for `scan_zip` when all structural zero indicators were equal to one.

# scanstatistics 0.1