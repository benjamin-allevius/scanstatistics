

# fast_kulldorff <- function(aggregates, 
#                            fk_priority_term_function,
#                            rel_tolerance = 0.01,
#                            max_iter = 100) {
#   foreach(W = unique(aggregates[, duration])) {
#     ags <- aggregates[duration == W]
#     # choose included streams at random
#     included_streams <- 1
#     # initialize relative risks at random
#     relative_risks <- 1
#     previous_score <- 1
#     for (i in seq(max_iter)) {
#       maxreg <- ags[stream %in% included_streams] %>%
#         fast_kulldorff_priority(relative_risks = ) %>%
#         fast_kulldorff_maxregion
#       
#       if (has_converged(maxreg[, score], previous_score, rel_tolerance)) {
#         maxreg
#       }
#       
#       relative_risks <- relative_risk_mle(ags, 
#                                           maxreg[, unlist(included_locations)])
#                 
#     }
#   }
# }

# G_W^D(s_i) for Poisson
fk_priority_fun_poisson <- function(c, b, q) {
  c * log(q) + b * (1 - q)
}

# G_W^D(s_i) for Gaussian
fk_priority_fun_gaussian <- function(c, b, q) {
  (q - 1) * (c - (q + 1) * b / 2)
}


# duration given, subset of streams given
# relative_risks is a data.table
fast_kulldorff_priority <- function(aggregates, 
                                    relative_risks,
                                    fk_priority_fun) {
  aggregates[, 
    .(priority = sum(fk_priority_fun(aggregate_count,
                                     aggregate_baseline,
                                     relative_risks[stream])),
                 included_streams = list(stream)),
             by = .(location, duration)]
}

#' Computes the region which maximizes the conditional expectation-based score.
#' 
#' This function computes the region \eqn{S*} of locations which maximize the
#' expectation-based score function, given the relative risks:
#' \deqn{
#' S* = argmax_S F(D,S,W|{q_m}) = argmax_S \sum_i G_W^D(s_i)
#' }
#' for a given set of streams \eqn{D} and a given event duration \eqn{W}.
#' Only those locations \eqn{s_i} which make a positive contribution to the sum,
#' i.e. those regions \eqn{s_i} for which \eqn{G_W^D(s_i) > 0}, are included.
#' @param priorities A \code{data.table} with columns \code{location, duration,
#'    priority, included_streams}. All elements in the colum \code{duration}
#'    should be the same, and the column \code{included_streams} is a list, in 
#'    which all elements are also equal.
#' @return A single-row \code{data.table} with columns \code{included_streams, 
#'    included_locations, duration, score}. The columns \code{included_streams}
#'    and \code{duration} have the same elements as in the input (now without
#'    duplicates). The column \code{included_locations} consists of those 
#'    locations which form the score-maximizing region. The column \code{score}
#'    contains the conditional expectation-based score; conditional on the 
#'    relative risks which entered into the calculations of the input 
#'    priorities.
fast_kulldorff_maxregion <- function(priorities) {
  priorities[priority > 0,
             .(included_streams = included_streams,
               included_locations = list(location),
               duration = duration, 
               score = sum(priority))][which.max(score)]
}

#' Compute the relative risk estimates for each stream by maximum likelihood.
#' 
#' This function computes the maximum likelihood relative risk estimates 
#' \eqn{q_m} for each data stream, as the ratio of aggregate counts 
#' \eqn{C^m(S,W)} to aggregate baselines \eqn{B^m(S,W)}, for a given region 
#' \eqn{S} and a given event duration \eqn{W}.
#' @param aggregates A \code{data.table} with columns \code{location, stream,
#'    duration, aggregate_count, aggregate_baseline}. All elements in the column
#'    \code{duration} should be the same.
#' @param locations A vector of locations, corresponding to the region \eqn{S}.
#' @return A vector of relative risks, one for each data stream. The vector is
#'    sorted in ascending order according to the data stream names (or numbers),
#'    and is named in that way.
relative_risk_mle <- function(aggregates, locations) {
  rr_mle <- aggregates[location %in% locations,
             .(aggregate_count = sum(aggregate_count),
               aggregate_baseline = sum(aggregate_baseline)),
             by = .(stream, duration)][, 
    .(relative_risk = max(1, aggregate_count / aggregate_baseline)),
    by = .(stream)]
  relrisks <- rr_mle[, relative_risk]
  names(relrisks) <- rr_mle[, stream]
  relrisks
}