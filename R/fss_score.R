
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
