

fast_kulldorff <- function(aggregates, 
                           fk_priority_term_function,
                           rel_tolerance = 0.01,
                           max_iter = 100) {
  foreach(W = unique(aggregates[, duration])) {
    ags <- aggregates[duration == W]
    # choose included streams at random
    included_streams <- 1
    # initialize relative risks at random
    relative_risks <- 1
    previous_score <- 1
    for (i in seq(max_iter)) {
      maxreg <- ags[stream %in% included_streams] %>%
        fast_kulldorff_priority(relative_risks = ) %>%
        fast_kulldorff_maxregion
      
      if (has_converged(maxreg[, score], previous_score, rel_tolerance)) {
        maxreg
      }
      
      relative_risks <- relative_risk_mle(ags, 
                                          maxreg[, unlist(included_locations)])
                
    }
  }
}


fk_priority_term_poisson <- function(c, b, q) {
  c * log(q) + b * (1 - q)
}

fk_priority_term_gaussian <- function(c, b, q) {
  (q - 1) * (c - (q + 1) * b / 2)
}


# duration given, subset of streams given
# relative_risks is a data.table
fast_kulldorff_priority <- function(aggregates, 
                                    relative_risks,
                                    fk_priority_term) {
  aggregates[, 
    .(priority = sum(fk_priority_term(aggregate_count,
                                      aggregate_baseline,
                                      relative_risks[stream])),
                 included_streams = list(stream)),
             by = .(location, duration)]
}

fast_kulldorff_maxregion <- function(fk_priorities) {
  fk_priorities[priority > 0,
                .(score = sum(priority),
                  included_streams = included_streams,
                  included_locations = list(location))][which.max(score)]
}

relative_risk_mle <- function(aggregates, locations) {
  rr_mle <- aggregates[location %in% locations,
             .(aggregate_count = sum(aggregate_count),
               aggregate_baseline = sum(aggregate_baseline)),
             by = .(stream, duration)][, 
    .(relative_risk = max(1, aggregate_count / aggregate_baseline)),
    by = .(stream)]
  rrs <- rr_mle[, relative_risk]
  names(rrs) <- rr_mle[, stream]
  rrs
}