# scanstatistics 1.0.2

## Minor changes

* Removed unneeded internal functions that caused package to not be loaded.
* Added `CITATION` file with updated citation for the package.
* Fixes to multiple bugs found by Kelly Reeve.

# scanstatistics 1.0.1

* Fixes to compiler warnings; needed to keep package on CRAN.

# scanstatistics 1.0

## Major changes

* New interface for main functions: accept data frames or matrices instead of
  data tables.
* All scan statistics reimplemented in C++.
* Several new scan statistics available.

### Minor changes

* The functions `knn_zones` and `flexible_zones` now run faster due to change
  in algorithms.

# scanstatistics 0.1
