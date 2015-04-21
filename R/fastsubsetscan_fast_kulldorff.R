
#' Determine the highest-scoring subsets of locations and streams by the Fast 
#' Kulldorff method.
#' 
#' Determine the highest-scoring region (subset of locations) and subset of data
#' streams by the Fast Kulldorff method.
#' @param counts A \code{data.table} with columns \code{time, location, stream,
#'    count, baseline}, and possibly other columns depending on the distribution
#'    used. For example, a column \code{variance} is needed for the normal 
#'    distribution. The column \code{count} is not limited to actual integer
#'    counts, but should instead be appropriate for the distribution used. The 
#'    column \code{time} could for example be POSIXct times, but must in any
#'    case be sortable so that the most recent time is the last element in the
#'    sorted structure when sorting into ascending order.
#' @param distribution One of "poisson", "gaussian" (or "normal"), or 
#'    "exponential".
#' @param restarts The number of random restarts to perform in finding the
#'    score-maximizing region. A higher number means a greater chance of finding
#'    a global rather than local maximum of the score function.
#' @param tol The tolerance used in \code{\link{has_converged}}.
#' @param max_iter The maximum number of iterations to perform in finding the
#'    score-maximizing region, if convergence of the score is not reached 
#'    sooner.
#' @return A \code{data.table} with columns \code{score, included_streams,
#'    region, duration}, keyed by \code{score} (so the last row contains the 
#'    highest score). If no positive-scoring subsets were found, the table will
#'    be empty.
fast_kulldorff <- function(counts, 
                           distribution = "poisson",
                           restarts = 50,
                           tol = 0.01,
                           max_iter = 100) {
  # Input preprocessing --------------------------------------------------------
  distribution <- tolower(distribution)
  restarts <- as.integer(restarts)
  max_iter <- as.integer(max_iter)
  
  # Input validation------------------------------------------------------------
  if (!is.data.table(counts)) {
    stop("The log-likelihoods must be supplied as a data.table.")
  }
  valid_distributions <- c("poisson", "gaussian", "normal", "exponential")
  if (distribution %notin% valid_distributions) {
    stop("distribution must be one of ",
         paste0(valid_distributions, collapse = ", "))
  }
  if (!all(c("stream", "location", "time", "count", "baseline") %in% 
             names(counts))) {
    stop("The data.table containing the counts must contain the columns ",
         "'stream', 'location', 'time', 'count', 'baseline'.")
  }
  if (distribution %in% c("normal", "gaussian") && 
        "variance" %notin% names(counts)) {
    stop("If distribution is 'normal' or 'gaussian', the counts table must ",
         "contain the column 'variance'.")
  }
  if (length(tol) != 1 || tol <= 0) {
    stop("tol must be a single number greater than 0.")
  }
  if (length(restarts) != 1 || length(max_iter) != 1 || 
        restarts < 1 || max_iter < 1) {
    stop("restarts and max_iter must be single positive integers.")
  }
  
  # Input valid; start calculations --------------------------------------------
  
  # Define aggregation and score functions based on distribution
  initial_aggregation_fun <- dispatch_initial_aggregation(distribution)
  score_fun <- dispatch_score_function(distribution)
  cond_score_fun <- dispatch_cond_score_function(distribution)
  
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
  setcolorder(res, c("score", "included_streams", "region", "duration"))
  # Return an empty table if no optimal regions and stream subsets were found.
  if (all(is.na(res[, score]))) {
    return(res[0])
  }
  unique(res[!is.na(score), .SD, keyby = .(score)])
}

#' Find the score-maximizing subsets for multiple initial values of the relative
#' risks.
#' 
#' This function finds the score-maximizing region and subset of data streams
#' for a given event duration, by multiple calls to 
#' \code{\link{find_maximizing_region}}. The initial values of the relative 
#' risks are chosen at random for each call, so that with multiple calls, the 
#' global maximum of the score function is more likely to have been found.
#' @param restarts The number of times we find the score-maximizing region and 
#'    subset of streams with different initial values.
#' @param ... Arguments passed to \code{\link{find_maximizing_subsets}}.
#' @return A \code{data.table} with \code{restarts} number of rows. This table
#'    has columns \code{duration, score, included_streams, region}.
random_restart_maximizer <- function(..., restarts = 50) {
  # Pre-allocate data.table
  maxed <- empty_optimal_stream_subset(nrow = restarts)
  for (i in seq(restarts)) {
    maxed[i, names(maxed) := find_maximizing_subsets(...)]
  }
  maxed
}

#' Find the score-maximizing subset of locations and data streams.
#' 
#' For a given event duration, find the subset of locations and data streams 
#' that maximize the given score function. 
#' @inheritParams relative_risk_mle
#' @inheritParams optimal_stream_subset
#' @param ... Arguments passed to \code{\link{find_maximizing_region}}.
#' @return A one-row \code{data.table} with columns \code{duration, score,
#'    included_streams, region}. If a score-maximizing region could not be 
#'    found, the table contains NAs. Otherwise, the colum \code{score} contains 
#'    the sum of the input scores over all data streams, for the given region 
#'    and duration, and the column \code{included_streams} contains those data 
#'    streams that made a positive contribution to this sum.
find_maximizing_subsets <- function(aggregates, score_fun, ...) {
  empty_output <- empty_optimal_stream_subset(nrow = 1)
  # Find the maximizing region by iterative procedure
  maxregion <- find_maximizing_region(aggregates, ...)
  if (length(maxregion) == 1 && is.na(maxregion)) {
    return(empty_output)
  }
  output <- optimal_stream_subset(aggregates = aggregates, 
                                  locations = maxregion,
                                  score_fun = score_fun)
  if (nrow(output) == 0) {
    return(empty_output)
  }
  output
}

#' Returns the minimal subset of streams that maximizes the score function.
#' 
#' For a given region (vector of locations) and event duration, return those 
#' data streams which make a positive contribution to the score, and return this 
#' score along with that minimal subset of streams that contribute to it.
#' @inheritParams relative_risk_mle
#' @param locations A vector of locations.
#' @param score_fun A two-parameter scalar input, single scalar output score 
#'    function.
#' @return A \code{data.table} with columns \code{duration, score, 
#'    included_streams, region}. The colum \code{score} contains the sum of the 
#'    input scores over all data streams, for the given region and duration. The 
#'    column \code{included_streams} contains those data streams that made a 
#'    positive contribution to this sum. The table has a single row if an 
#'    optimal subset of streams was found; it has zero rows (is empty) if an
#'    optimal subset of streams was not found.
optimal_stream_subset <- function(aggregates, locations, score_fun) {
  aggregates %>%
    aggregate_per_stream(locations = locations, TRUE) %>%
    expectation_based_score(score_function = score_fun, TRUE) %>%
    score_minimal_stream_subset(region_as_list = TRUE)
}

# Performs the two-step iterative procedure to find conditionally optimal region
# Generates initial random relative risks
# Outputs a vector of locations
#' Find the conditional score-maximizing region.
#' 
#' Find the region that maximizes the conditional score function by iteratively
#' finding the subset of locations (the region) that maximizes the score given
#' the current relative risk, and finding the maximum likelihood estimates of 
#' the relative risks given the region.
#' @inheritParams relative_risk_mle
#' @inheritParams fast_kulldorff_priority
#' @param max_iter The maximum number of iterations before returning the current
#'    best candidate for the score-maximizing region. 
#' @param ... Parameters sent to \code{\link{has_converged}}.
#' @return A vector of locations, or NA if no such optimal region is found. The
#'    function returns the best candidate for the score-maximizing region if 
#'    convergence of the score is reached before \code{max_iter} iterations. If
#'    this convergence is not obtained before \code{max_iter} iterations, then
#'    the current best candidate (which may be NA) for the score-maximizing 
#'    region is returned.
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
                              cond_score_fun = cond_score_fun)
    
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
#' included by this, choose one stream uniformly at random to return.
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
#' @param cond_score_fun The score function, conditional on known relative 
#'    risks. One for each data stream; corresponds to a term in the sum of the
#'    priority \eqn{G_W^D(s_i)}. Takes three scalar inputs: an aggregate count, 
#'    an aggregate baseline, and a relative risk. Outputs a scalar value, the 
#'    conditional score.
#' @return A \code{data.table} with columns \code{location, duration, priority,
#'    included_streams}.
fast_kulldorff_priority <- function(aggregates, 
                                    relative_risks,
                                    cond_score_fun) {
  aggregates[, 
    .(included_streams = list(stream),
      priority = sum(cond_score_fun(aggregate_count,
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