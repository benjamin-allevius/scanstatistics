
#' Calculate the aggregate counts and baselines for the EBP scan statistic
#' over all event durations.
#' 
#' Calculate the aggregate counts \eqn{C_{i,m}} and aggregate baselines 
#' \eqn{B_{i,m}} for the expectation-based Poisson (EBP) scan statistic, for 
#' each location \eqn{i} and data stream \eqn{m}. I.e. the cumulative sum over
#' the event duration, from shortest to longest, is calculated.
#' 
#' @param counts A \code{data.table} with columns \code{location, stream, 
#'        duration, count, baseline}.
aggregate_CB_poisson <- function(counts) {
  counts[, .(duration = duration,
             aggregate_count = cumsum(count),
             aggregate_baseline = cumsum(baseline)),
         by = .(location, stream)]
}

#' Calculate the aggregate counts and baselines for the EBG scan statistic
#' over all event durations.
#' 
#' Calculate the aggregate counts \eqn{C_{i,m}} and aggregate baselines 
#' \eqn{B_{i,m}} for the expectation-based Gaussian (EBG) scan statistic, for 
#' each location \eqn{i} and data stream \eqn{m}. I.e. the cumulative sum over
#' the event duration, from shortest to longest, is calculated.
#' 
#' @inheritParams aggregate_CB_poisson
aggregate_CB_gaussian <- function(counts) {
  counts[, .(duration = duration,
             aggregate_count = cumsum(count * baseline / variance),
             aggregate_baseline = cumsum(baseline^2 / variance)),
         by = .(location, stream)]
}

#' Calculate the aggregate counts and baselines for the EBE scan statistic
#' over all event durations.
#' 
#' Calculate the aggregate counts \eqn{C_{i,m}} and aggregate baselines 
#' \eqn{B_{i,m}} for the expectation-based Exponential (EBE) scan statistic, for 
#' each location \eqn{i} and data stream \eqn{m}. I.e. the cumulative sum over
#' the event duration, from shortest to longest, is calculated.
#' 
#' @inheritParams aggregate_CB_poisson
aggregate_CB_exponential <- function(counts) {
  counts[, .(duration = duration,
             aggregate_count = cumsum(count / baseline),
             aggregate_baseline = duration), 
         by = .(location, stream)]
}

#' Calculate the count and baseline aggregates over each region-stream-duration
#' combination.
#' 
#' Take the already calculated aggregates for each location, stream and duration
#' and sum them over all locations in each region.
#' 
#' @param aggregates A \code{data.table} containing columns \code{region, 
#'        duration, stream, location, aggregate_count, aggregate_baseline}.
#'        The latter two columns contain the aggregate counts and baselines
#'        as produced by the functions \code{aggregate_CB_X}, where \code{X} is
#'        e.g. \code{poisson} or \code{gaussian}.
#' @return A \code{data.table} with the same columns except \code{location}; the
#'         aggregate quantities have now been summed over all locations in each
#'         region, for each region, duration, and data stream.
aggregate_again <- function(aggregates) {
  aggregates[, .(aggregate_count = sum(aggregate_count),
                 aggregate_baseline = sum(aggregate_baseline)),
             by = .(region, duration, stream)]
}


#' Calculates the score for the expectation-based Poisson scan statistic.
#' 
#' Calculates the score for the expectation-based Poisson scan statistic, given
#' scalar aggregate counts and baselines.
#' 
#' @param c A scalar; an aggregate count.
#' @param b A scalar; an aggregate baseline.
score_fun_EBP <- function(c, b) {
  ifelse(c > b, c * (log(c) - log(b) - 1) + b, 0)
}

#' Calculates the expectation-based Poisson score for each region, stream, and
#' duration combination.
#' 
#' Given the aggegate counts and baselines for each region-stream-duration 
#' combination, calculates the expectation-based Poisson score function value.
#' 
#' @param aggregates A \code{data.table} with columns \code{region, duration, 
#'        stream, aggregate_count, aggregate_baseline}.
#' @return A \code{data.table} with columns \code{region, duration, 
#'         stream, score}.
score_EBP <- function(aggregates) {
  aggregates[, .(score = score_fun_EBP(aggregate_count, aggregate_baseline)),
             by = .(region, duration, stream)]
}

#' Score and minimal stream subset for the Naive Kulldorff method.
#' 
#' Calculates the score for each combination of region and duration, and the 
#' minimal subset of locations that contribute to the score, according to the
#' Naive Kulldorff method. See Neill et. al. (2013) section 3.2.
#' 
#' @scores A \code{data.table} with columns \code{region, duration, stream, 
#'         score}.
#' @return A \code{data.table} with columns \code{region, duration, stream, 
#'         score, included_streams}. The colum \code{score} contain the sum
#'         of the input score over all data streams, for each region and 
#'         duration. The column \code{included_streams} contain those data 
#'         streams that made a positive contribution to this sum.
score_minimal_stream_subset <- function(scores) {
  scores[score > 0, 
         .(score = sum(score), included_streams = list(stream)), 
         by = .(region, duration)]
}


#' Calculate the scores according to the Naive Kulldorff method.
#' @importFrom magrittr %>%
naive_kulldorff_poisson <- function(counts, regions) {
  counts %>%
    aggregate_CB_poisson %>%
    region_joiner(regions = regions, 
                  keys = c("region", "duration", "stream")) %>%
    aggregate_again %>%
    score_EBP %>%
    score_minimal_stream_subset
}
