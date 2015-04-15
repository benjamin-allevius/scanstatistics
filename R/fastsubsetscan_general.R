# General functions ------------------------------------------------------------

#' Calculate the count and baseline aggregates over each region-stream-duration
#' combination.
#' 
#' Take the already calculated aggregates for each location, stream and duration
#' and sum them over all locations in each region. These are the quantities 
#' denoted \eqn{C^m(S,W)} and \eqn{B^m(S,W)}.
#' @param aggregates See \code{\link{aggregate_per_stream_regionnotlist}} if
#'    argument \code{region_as_list} is \code{TRUE}. See 
#'    \code{\link{aggregate_per_stream_regionaslist}} if argument 
#'    \code{region_as_list} is \code{FALSE}.
#' @param location A vector of locations corresponding to a single region; 
#'    passed to \code{\link{aggregate_per_stream_regionaslist}} if 
#'    \code{region_as_list} is \code{TRUE}.
#' @param region_as_list Boolean: if \code{TRUE}, call 
#'    \code{\link{aggregate_per_stream_regionaslist}} with input. If 
#'    \code{FALSE}, call \code{\link{aggregate_per_stream_regionnotlist}} with 
#'    input.
aggregate_per_stream <- function(aggregates, 
                                 locations = NULL, 
                                 region_as_list = FALSE) {
  if (region_as_list) {
    return(aggregate_per_stream_regionaslist(aggregates, locations))
  } else {
    return(aggregate_per_stream_regionnotlist(aggregates))
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
aggregate_per_stream_regionaslist <- function(aggregates, locations) {
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
aggregate_per_stream_regionnotlist <- function(aggregates) {
  aggregates[, .(aggregate_count = sum(aggregate_count),
                 aggregate_baseline = sum(aggregate_baseline)),
             by = .(region, duration, stream)]
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
score_EB <- function(aggregates, score_function) {
  aggregates[, .(score = score_function(aggregate_count, aggregate_baseline)),
             by = .(region, duration, stream)]
}

# Expectation-based Poisson ----------------------------------------------------

#' Calculate the aggregate counts and baselines for the EBP scan statistic
#' over all event durations.
#' 
#' Calculate the aggregate counts \eqn{C_{i,m}(W)} and aggregate baselines 
#' \eqn{B_{i,m}(W)} for the expectation-based Poisson (EBP) scan statistic, for 
#' each location \eqn{i}, data stream \eqn{m}, and event duration \eqn{W}. In 
#' essence, the cumulative sum over the event duration, from shortest to 
#' longest, is calculated.
#' @param counts A \code{data.table} with columns \code{location, stream, 
#'    duration, count, baseline}.
aggregate_CB_poisson <- function(counts) {
  counts[, .(duration = duration,
             aggregate_count = cumsum(count),
             aggregate_baseline = cumsum(baseline)),
         by = .(location, stream)]
}

#' Calculates the score for the expectation-based Poisson scan statistic.
#' 
#' Calculates the score for the expectation-based Poisson scan statistic, given
#' scalar aggregate counts and baselines.
#' @param c A scalar; an aggregate count.
#' @param b A scalar; an aggregate baseline.
score_fun_EBP <- function(c, b) {
  ifelse(c > b, c * (log(c / b) - 1) + b, 0)
}

# Expectation-based Gaussian ---------------------------------------------------

#' Calculate the aggregate counts and baselines for the EBG scan statistic
#' over all event durations.
#' 
#' Calculate the aggregate counts \eqn{C_{i,m}} and aggregate baselines 
#' \eqn{B_{i,m}} for the expectation-based Gaussian (EBG) scan statistic, for 
#' each location \eqn{i} and data stream \eqn{m}. I.e. the cumulative sum over
#' the event duration, from shortest to longest, is calculated.
#' @inheritParams aggregate_CB_poisson
aggregate_CB_gaussian <- function(counts) {
  counts[, .(duration = duration,
             aggregate_count = cumsum(count * baseline / variance),
             aggregate_baseline = cumsum(baseline^2 / variance)),
         by = .(location, stream)]
}

#' Calculates the score for the expectation-based Gaussian scan statistic.
#' 
#' Calculates the score for the expectation-based Gaussian scan statistic, given
#' scalar aggregate counts and baselines.
#' @param c A scalar; an aggregate count.
#' @param b A scalar; an aggregate baseline.
score_fun_EBG <- function(c, b) {
  ifelse(c > b, (c - b)^2 / (2 * b), 0)
}

# Expectation-based Exponential ------------------------------------------------

#' Calculate the aggregate counts and baselines for the EBE scan statistic
#' over all event durations.
#' 
#' Calculate the aggregate counts \eqn{C_{i,m}} and aggregate baselines 
#' \eqn{B_{i,m}} for the expectation-based Exponential (EBE) scan statistic, for 
#' each location \eqn{i} and data stream \eqn{m}. I.e. the cumulative sum over
#' the event duration, from shortest to longest, is calculated.
#' @inheritParams aggregate_CB_poisson
aggregate_CB_exponential <- function(counts) {
  counts[, .(duration = duration,
             aggregate_count = cumsum(count / baseline),
             aggregate_baseline = duration), 
         by = .(location, stream)]
}

#' Calculates the score for the expectation-based exponential scan statistic.
#' 
#' Calculates the score for the expectation-based exponential scan statistic, 
#' given scalar aggregate counts and baselines.
#' @param c A scalar; an aggregate count.
#' @param b A scalar; an aggregate baseline.
score_fun_EBE <- function(c, b) {
  ifelse(c > b, b * (log(b / c) - 1) + c, 0)
}
