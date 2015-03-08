





logsumexp_llh_over_regions <- function(spatial_llhs) {
  spatial_llhs[, .(llh = logsumexp(llh)), by = .(event)]
}