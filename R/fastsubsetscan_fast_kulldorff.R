

fast_kulldorff <- function(counts, 
                           distribution = "poisson",
                           restarts = 50,
                           tol = 0.01,
                           max_iter = 100) {
  # [Input validation here]
  
  # Define aggregation and score functions based on distribution
  initial_aggregation_fun <- aggregate_CB(distribution)
  score_fun <- score_function_EB(distribution)
  cond_score_fun <- conditional_score_function_EB(distribution)
  
  # Do initial aggregation for each stream, location, and duration
  aggregates <- counts %>%
    add_duration(keys = c("duration", "location", "stream")) %>%
    initial_aggregation_fun
  
  # For each duration, find score-maximizing subsets of locations and streams
  durations <- unique(aggregates[, duration])
  res <- foreach::foreach(W = durations, 
                          .combine = rbind, 
                          .inorder = FALSE) %do% {
    ags <- aggregates[duration == W]
    random_restart_maximizer(aggregates = ags,
                             score_fun = score_fun,
                             cond_score_fun = cond_score_fun,
                             tol = tol,
                             max_iter = max_iter,
                             restarts = restarts)
  }
  unique(res[!is.na(duration), .SD, keyby = .(score)])
}

random_restart_maximizer <- function(..., restarts = 50) {
  # Pre-allocate data.table
  maxed <- empty_max_table(nrow = restarts)
  for (i in seq(restarts)) {
    maxed[i, names(maxed) := find_maximizing_subsets(...)]
  }
  maxed
}

# Table as returned by optimal_stream_subset
empty_max_table <- function(nrow = 1) {
  data.table(duration = rep(as.integer(NA), nrow), 
             score = as.numeric(NA), 
             included_streams = list(), 
             region = list())
}

# called for each duration W
# i.e. pass in aggregates for a single W
find_maximizing_subsets <- function(aggregates, score_fun, ...) {  
  # Find the maximizing region by iterative procedure
  maxregion <- find_maximizing_region(aggregates, ...)
  if (length(maxregion) == 1 && is.na(maxregion)) {
    return(empty_max_table(nrow = 1))
  }
  optimal_stream_subset(aggregates = aggregates, 
                        locations = maxregion,
                        score_fun = score_fun)
}

optimal_stream_subset <- function(aggregates, locations, score_fun) {
  aggregates %>%
    aggregate_per_stream(locations = locations, TRUE) %>%
    expectation_based_score(score_function = score_fun, TRUE) %>%
    score_minimal_stream_subset(region_as_list = TRUE)
}

# Performs the two-step iterative procedure to find conditionally optimal region
# Generates initial random relative risks
# Outputs a vector of locations
find_maximizing_region <- function(aggregates, 
                                   cond_score_fun, 
                                   max_iter = 100,
                                   ...) {
  # Choose which streams to initially include in maximization
  all_streams <- get_all_streams(aggregates)
  incl_streams <- choose_streams_randomly(all_streams)
  
  # Randomly initialize relative risks for included streams
  rel_risks <- random_relative_risk(incl_streams, all_streams)
   
  # Iterate until score is maximized
  current_maxregion <- as.integer(NA)
  previous_score <- as.numeric(NA)
  for (i in seq(max_iter)) {
    # For the current relative risks, find the region which maximizes score
    maxregscore <- aggregates[stream %in% incl_streams] %>%
      fast_kulldorff_priority(relative_risks = rel_risks,
                              conditional_score = cond_score_fun)
    
    # If all priorities are negative, need to reset relative risks and streams
    if (all(maxregscore[, priority] < 0)) {
      incl_streams <- choose_streams_randomly(all_streams)
      rel_risks <- random_relative_risk(incl_streams, all_streams)
      previous_score <- as.numeric(NA)
      next
    }
    maxregscore <- fast_kulldorff_maxregion(maxregscore)
    
    # Extract locations in region as vector
    current_maxregion <- maxregscore[, region][[1]]
    
    if (!is.na(previous_score) && 
          has_converged(maxregscore[, score], previous_score, ...)) {
      return(current_maxregion)
    }
    rel_risks <- relative_risk_mle(aggregates, current_maxregion)
    incl_streams <- all_streams[rel_risks > 1]
    previous_score <- maxregscore[, score]
  }
  # Did not converge before maximum number of iterations, but return anyway
  current_maxregion
}

#' Randomly assign relative risks for included streams, and set to 1 for 
#' non-included streams.
#' 
#' For the included streams, set the relative risks to a random number greater
#' than one. For the non-included streams, set the relative risk to one.
#' @param incl_streams A subset of the atomic vector all_streams.
#' @param all_streams An atomic vector containing the sorted names or 
#'    identifiers for all data streams.
random_relative_risk <- function(incl_streams, all_streams) {
  n_streams <- length(all_streams)
  ifelse(all_streams %in% incl_streams, 
         exp(2 * runif(n_streams)), 
         rep(1, n_streams))
}


#' Decide which data streams are to be included by independent coin flips.
#' 
#' Given a vector of data stream names, return a subset of these by including
#' each data stream based on the flip of a biased coin; the bias is random and 
#' the coin flips are independent. However, if by chance no streams are 
#' included, choose one stream uniformly at random to return.
#' @param all_streams A vector of data stream names or other identifiers.
choose_streams_randomly <- function(all_streams) {
  n_streams <- length(all_streams)
  
  # Decide which streams to include by coin flip
  stream_included <- rbinom(n_streams, 1, runif(1)) == 1
  
  # Make sure at least 1 stream is included
  if (all(!stream_included)) {
    stream_included[sample(seq(stream_included), 1)] <- TRUE
  }
  all_streams[stream_included]
}

#' Calculates the priority for each location, for a given duration and subset of
#' streams.
#' 
#' For a given duration and set of data streams, this function calculates the 
#' priority \eqn{G_W^D(s_i)} for each location, using the given relative risks
#' and conditonal score function (conditional on relative risks).
#' @inheritParams relative_risk_mle
#' @param relative_risks A vector of relative risks, one for each data stream.
#' @param conditional_score The score function, conditional on known relative 
#'    risks. One for each data stream; corresponds to a term in the sum of the
#'    priority \eqn{G_W^D(s_i)}. Takes three scalar inputs: an aggregate count, 
#'    an aggregate baseline, and a relative risk. Outputs a scalar value, the 
#'    conditional score.
#' @return A \code{data.table} with columns \code{location, duration, priority,
#'    included_streams}.
fast_kulldorff_priority <- function(aggregates, 
                                    relative_risks,
                                    conditional_score) {
  aggregates[, 
    .(included_streams = list(stream),
      priority = sum(conditional_score(aggregate_count,
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