# Note about potential functionality duplication: the data.table package does 
# not allow columns that are lists to be used as keys, which means you can't
# do 'by' groupings on these columns. Therefore, we need one function for the
# case in which the column is of atomic type, and another for when the column is 
# a list.
# Specific consequence: need to provide separate functions for the case in which
# a region is given in a data.table as an atomic vector (integer most likely)
# and for the case in which the region column is given as a list, each item a
# vector of locations.


# score_minimal_stream_subset --------------------------------------------------

#' Score and minimal stream subset for the Kulldorff method.
#' 
#' Calculate the score and the subset of data streams that make a (positive)
#' contribution to this score. Calls either 
#' \code{\link{score_minimal_stream_subset_list}} or
#' \code{\link{score_minimal_stream_subset_atomic}}.
#' @param scores A \code{data.table} containing expectation-based scores, which
#'    are to be summed over all streams with positive scores.
#' @inheritParams aggregate_per_stream
score_minimal_stream_subset <- function(scores, region_as_list = FALSE) {
  if (region_as_list) {
    return(score_minimal_stream_subset_list(scores))
  } else {
    return(score_minimal_stream_subset_atomic(scores))
  }
}

#' Score and minimal stream subset for the Kulldorff method.
#' 
#' Calculates the score for a given region and event duration, by summing those
#' per-stream scores which are positive. These streams are returned in the 
#' output.
#' @param scores A \code{data.table} with columns \code{region, duration, 
#'    stream, score}. The column \code{region} should be a list with all 
#'    elements equal, each element the same vector of locations that constitute 
#'    the region. Likewise, the column \code{duration} should have all elements
#'    equal.
#' @return A \code{data.table} with columns \code{region, duration, 
#'    included_streams, score}. The colum \code{score} contain the sum of the 
#'    input scores over all data streams, for each region and duration. The 
#'    column \code{included_streams} contain those data streams that made a 
#'    positive contribution to this sum.
score_minimal_stream_subset_list <- function(scores) {
  reg <- scores[1, region]
  scores[score > 0, 
         .(score = sum(score), included_streams = list(stream)), 
         by = .(duration)][, region := list(reg)]
}

#' Score and minimal stream subset for the Kulldorff method.
#' 
#' Calculates the score for each combination of region and duration, and the 
#' minimal subset of locations that contribute to the score, according to the
#' Kulldorff method. See Neill et. al. (2013) section 3.2.
#' @param scores A \code{data.table} with columns \code{region, duration, stream, 
#'    score}.
#' @return A \code{data.table} with columns \code{region, duration, stream, 
#'    score, included_streams}. The colum \code{score} contain the sum of the 
#'    input score over all data streams, for each region and duration. The 
#'    column \code{included_streams} contain those data streams that made a 
#'    positive contribution to this sum.
score_minimal_stream_subset_atomic <- function(scores) {
  scores[score > 0, 
         .(included_streams = list(stream), score = sum(score)), 
         by = .(region, duration)]
}

# aggregate_per_stream ---------------------------------------------------------


#' Calculate the count and baseline aggregates over each region-stream-duration
#' combination.
#' 
#' Take the already calculated aggregates for each location, stream and duration
#' and sum them over all locations in each region. These are the quantities 
#' denoted \eqn{C^m(S,W)} and \eqn{B^m(S,W)}.
#' @param aggregates See \code{\link{aggregate_per_stream_regionnotlist}} if
#'    argument \code{region_as_list} is \code{TRUE}. See 
#'    \code{\link{aggregate_per_stream_list}} if argument 
#'    \code{region_as_list} is \code{FALSE}.
#' @param location A vector of locations corresponding to a single region; 
#'    passed to \code{\link{aggregate_per_stream_list}} if 
#'    \code{region_as_list} is \code{TRUE}.
#' @param region_as_list Boolean: is the region colum a list, or an atomic 
#'    vector? 
aggregate_per_stream <- function(aggregates, 
                                 locations = NULL, 
                                 region_as_list = FALSE) {
  if (region_as_list) {
    return(aggregate_per_stream_list(aggregates, locations))
  } else {
    return(aggregate_per_stream_atomic(aggregates))
  }
}

#' Compute the per-stream aggregates for a given region and event duration.
#' 
#' Compute the per-stream aggregate counts \eqn{C^m(S,W)} and aggregate 
#' baselines \eqn{B^m(S,W)}, given a single region (specified as a set of 
#' locations) and a single duration.
#' @inheritParams relative_risk_mle
#' @return A new \code{data.table} with columns \code{region, stream, duration,
#'    aggregate_count, aggregate_baseline}. The column \code{region} is a list
#'    which contains identical elements; the column duration also contains 
#'    identical elements.
aggregate_per_stream_list <- function(aggregates, locations) {
  aggregates[location %in% locations,
             .(region = list(location),
               aggregate_count = sum(aggregate_count),
               aggregate_baseline = sum(aggregate_baseline)),
             by = .(stream, duration)]
}

#' Calculate the count and baseline aggregates over each region-stream-duration
#' combination.
#' 
#' Take the already calculated aggregates for each location, stream and duration
#' and sum them over all locations in each region. These are the quantities 
#' denoted \eqn{C^m(S,W)} and \eqn{B^m(S,W)}.
#' @param aggregates A \code{data.table} containing columns \code{region, 
#'    duration, stream, location, aggregate_count, aggregate_baseline}. The 
#'    latter two columns contain the aggregate counts and baselines as produced 
#'    by the functions \code{aggregate_CB_X}, where \code{X} is e.g. 
#'    \code{poisson} or \code{gaussian}.
#' @return A \code{data.table} with the same columns except \code{location}; the
#'    aggregate quantities have now been summed over all locations in each 
#'    region, for each region, duration, and data stream.
aggregate_per_stream_atomic <- function(aggregates) {
  aggregates[, .(aggregate_count = sum(aggregate_count),
                 aggregate_baseline = sum(aggregate_baseline)),
             by = .(region, duration, stream)]
}

# Expectation-based score for data.table (not scalar input) =-------------------

expectation_based_score <- function(aggregates, 
                                    score_function, 
                                    region_as_list = FALSE) {
  if (region_as_list) {
    return(expectation_based_score_list(aggregates, score_function))
  } else {
    return(expectation_based_score_atomic(aggregates, score_function))
  }  
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
expectation_based_score_list <- function(aggregates, score_function) {
  aggregates[, .(region = region,
                 score = score_function(aggregate_count, aggregate_baseline)),
             by = .(duration, stream)]
}

#' Calculates the expectation-based score for each region, stream, and duration 
#' combination.
#' 
#' Given the aggegate counts and baselines for each region-stream-duration 
#' combination, calculates the expectation-based score function value, with the
#' given score function.
#' @param aggregates A \code{data.table} with columns \code{region, duration, 
#'    stream, aggregate_count, aggregate_baseline}.
#' @param score_function A two-parameter scalar input, single scalar output
#'    score function.
#' @return A \code{data.table} with columns \code{region, duration, stream, 
#'    score}.
expectation_based_score_atomic <- function(aggregates, score_function) {
  aggregates[, .(score = score_function(aggregate_count, aggregate_baseline)),
             by = .(region, duration, stream)]
}
