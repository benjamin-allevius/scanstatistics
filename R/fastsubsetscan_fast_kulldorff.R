

fast_kulldorff <- function(counts, 
                           distribution = "poisson",
                           random_restarts = 50,
                           tol = 0.01,
                           max_iter = 100) {
  # Input validation here
  initial_aggregation_fun <- aggregate_CB(distribution)
  score_fun <- score_function_EB(distribution)
  conditional_score_fun <- conditional_score_function_EB(distribution)
  
  aggregates <- counts %>%
    add_duration(keys = c("duration", "location", "stream")) %>%
    initial_aggregation_fun
  # Choose priority and score functions based on distribution
  durations <- unique(aggregates[, duration])
  foreach::foreach(W = durations, .combine = rbind, .inorder = FALSE) %do% {
    ags <- aggregates[duration == W]
    random_restart_maximizer(ags)
  }
}

random_restart_maximizer <- function(aggregates,
                                     priority_term, 
                                     score_fun,
                                     restarts = 50,
                                     ...) {
  foreach::foreach(i = seq(restarts),
                  .combine = rbind,
                  .inorder = FALSE) %do% {
    find_maximizing_subsets(aggregates, priority_term, score_fun, ...)
  }
}

# called for each duration W
# i.e. pass in aggregates for a single W
find_maximizing_subsets <- function(aggregates, 
                                    priority_term, 
                                    score_fun,
                                    ...) {
  all_streams <- sort(unique(aggregates[, stream]))
  n_streams <- length(all_streams)
  # Decide which streams to include by coin flip
  stream_included <- rbinom(n_streams, 1, runif(1)) == 1
  # Make sure at least 1 stream is included
  if (all(!stream_included)) {
    stream_included[sample(seq(stream_included), 1)] <- TRUE
  }
  incl_streams <- all_streams[stream_included]
  # Randomly initialize relative risks for included streams
  rel_risks <- ifelse(stream_included, 
                      exp(2 * runif(n_streams)), 
                      rep(1, n_streams))
  # Iterative maximization step
  maxregion <- find_maximizing_region(aggregates, 
                                      all_streams, 
                                      incl_streams, 
                                      priority_term,
                                      rel_risks, 
                                      tol, 
                                      max_iter)
  # Region is conditionally optimal given relative risks
  maxreg_locations <- maxregion[, region][[1]]
  # Find the optimal subset of streams for the region found above
  optimal_subset <- aggregates %>%
    aggregate_per_stream(locations = maxreg_locations,
                         region_as_list = TRUE) %>%
    expectation_based_score(score_function = score_fun,
                            region_as_list = TRUE) %>%
    score_minimal_stream_subset(region_as_list = TRUE)
}

# Performs the two-step iterative procedure to find conditionally optimal region
find_maximizing_region <- function(aggregates, 
                                   all_streams,
                                   incl_streams, 
                                   priority_term,
                                   rel_risks, 
                                   tol, 
                                   max_iter) {
  previous_score <- 1
  for (i in seq(max_iter)) {
    # Iterate until score is maximized
    maxregion <- aggregates[stream %in% incl_streams] %>%
      fast_kulldorff_priority(relative_risks = rel_risks,
                              priority_term = priority_term) %>%
      fast_kulldorff_maxregion
    
    rel_risks <- aggregates %>% 
      relative_risk_mle(locations = maxregion[, region][[1]])
    
    if (has_converged(maxregion[, score], previous_score, rel_tolerance)) {
      return(maxregion)
    }
    incl_streams <- all_streams[rel_risks > 1]
    previous_score <- maxregion[, score]
  }
  # Did not converge before maximum number of iterations, but return anyway
  maxregion
}

#' Calculates the priority for each location, for a given duration and subset of
#' streams.
#' 
#' For a given duration and set of data streams, this function calculates the 
#' priority \eqn{G_W^D(s_i)} for each location, using the input relative risks
#' and priority function.
#' @inheritParams relative_risk_mle
#' @param relative_risks A vector of relative risks, one for each data stream.
#' @param priority_term A function for the term in the sum over streams for the
#'    priority function \eqn{G_W^D(s_i)}. This function should take three 
#'    arguments: an aggregate count, an aggregate baseline, and a relative risk.
#' @return A \code{data.table} with columns \code{location, duration, priority,
#'    included_streams}.
fast_kulldorff_priority <- function(aggregates, 
                                    relative_risks,
                                    priority_term) {
  aggregates[, 
    .(included_streams = list(stream),
      priority = sum(priority_term(aggregate_count,
                                   aggregate_baseline,
                                   relative_risks[stream]))),
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
#'    region, duration, score}. The columns \code{included_streams}
#'    and \code{duration} have the same elements as in the input (now without
#'    duplicates). The column \code{region} consists of those 
#'    locations which form the score-maximizing region. The column \code{score}
#'    contains the conditional expectation-based score; conditional on the 
#'    relative risks which entered into the calculations of the input 
#'    priorities.
fast_kulldorff_maxregion <- function(priorities) {
  priorities[priority > 0,
             .(included_streams = included_streams,
               region = list(location),
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
    by = .(stream, duration)]
  relrisks <- rr_mle[, relative_risk]
  names(relrisks) <- rr_mle[, stream]
  relrisks
}