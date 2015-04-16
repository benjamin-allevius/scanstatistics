

fast_kulldorff <- function(aggregates, 
                           distribution = "poisson",
                           random_restarts = 50,
                           tol = 0.01,
                           max_iter = 100) {
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
    single_region_score_EB(score_function = score_fun) %>%
    single_region_minimal_stream_subset
    
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

#' Calculates the expectation-based score for each stream and the given region
#' and event duration.
#' 
#' Given the aggegate counts and baselines for each data stream and the given
#' region and event duration, calculates the expectation-based score, with the
#' given score function.
#' @param aggregates A \code{data.table} with columns \code{region, duration, 
#'    stream, aggregate_count, aggregate_baseline}.
#' @param score_function A two-parameter scalar input, single scalar output
#'    score function.
#' @return A \code{data.table} with columns \code{region, duration, stream, 
#'    score}.
single_region_score_EB <- function(aggregates, score_function) {
  aggregates[, .(region = region,
                 score = score_function(aggregate_count, aggregate_baseline)),
             by = .(duration, stream)]
}

#' Score and minimal stream subset for the Kulldorff method.
#' 
#' Calculates the score for each combination of region and duration, and the 
#' minimal subset of locations that contribute to the score, according to the
#' Kulldorff method. See Neill et. al. (2013) section 3.2.
#' @param scores A \code{data.table} with columns \code{region, duration, 
#'    stream, score}. The column \code{region} should be a list with all 
#'    elements equal, each element the same vector of locations that constitute 
#'    the region. Likewise, the column \code{duration} should have all elements
#'    equal.
#' @return A \code{data.table} with columns \code{region, duration, stream, 
#'    score, included_streams}. The colum \code{score} contain the sum of the 
#'    input scores over all data streams, for each region and duration. The 
#'    column \code{included_streams} contain those data streams that made a 
#'    positive contribution to this sum.
single_region_minimal_stream_subset <- function(scores) {
  # Can't sum by list column in data.table, so must do workaround.
  reg <- scores[1, region]
  topscore <- scores[score > 0, 
                     .(score = sum(score), 
                       included_streams = list(stream)), 
                     by = .(duration)][, region := reg]
  setcolorder(topscore, c("included_streams", 
                          "region", 
                          "duration",
                          "score"))
  topscore
}


#' Is the relative error between two numbers is less than the given tolerance?
#' 
#' Given two consecutive numbers in a sequence, return \code{TRUE} if the
#' relative change is positive but less than the given tolerance.
#' @param current A scalar; the most recent value of the sequence.
#' @param previous A scalar; the second most recent value of the sequence, or a
#'    reference value.
#' @param tol The tolerance, a positive scalar near zero.
has_converged <- function(current, previous, tol = 0.01) {
  rel_change <- (current - previous) / abs(previous)
  rel_change > 0 && rel_change < tol
}

# term for a given stream in the G_W^D(s_i) sum for Poisson
priority_term_poisson <- function(c, b, q) {
  c * log(q) + b * (1 - q)
}

# term for a given stream in the G_W^D(s_i) sum for Gaussian
priority_term_gaussian <- function(c, b, q) {
  (q - 1) * (c - (q + 1) * b / 2)
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