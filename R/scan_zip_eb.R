

scan_zip_eb <- function(counts,
                        zones,
                        baselines = NULL,
                        probs = NULL,
                        n_mcsim = 0,
                        gumbel = TRUE, 
                        rel_tol = 1e-2) {
  
  zones_flat <- unlist(zones)
  zone_lengths <- unlist(lapply(zones, length))
  
  
}