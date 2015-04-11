
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