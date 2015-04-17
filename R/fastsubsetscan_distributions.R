
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

#' Expectation-based Poisson score, conditional on the relative risk.
#' 
#' This function computes the expectation-based score for Poisson-distributed 
#' counts, conditional on the relative risk. Appears as a term (corresponding to
#' a single data stream) in the sum for the priority function G_W^D(s_i) used in
#' the Fast Kulldorff method.
#' @param c Scalar; an aggregate count.
#' @param b Scalar; an aggregate baseline.
#' @param q Scalar; a relative risk.
conditional_score_fun_EBP <- function(c, b, q) {
  c * log(q) + b * (1 - q)
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

#' Expectation-based Gaussian score, conditional on the relative risk.
#' 
#' This function computes the expectation-based score for normal-distributed 
#' counts, conditional on the relative risk. Appears as a term (corresponding to
#' a single data stream) in the sum for the priority function G_W^D(s_i) used in
#' the Fast Kulldorff method.
#' @inheritParams conditional_score_fun_EBP
conditional_score_fun_EBG <- function(c, b, q) {
  (q - 1) * (c - (q + 1) * b / 2)
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

#' Expectation-based exponential-distribution score, conditional on the relative 
#' risk.
#' 
#' This function computes the expectation-based score for exponentially 
#' distributed counts, conditional on the relative risk. Appears as a term 
#' (corresponding to a single data stream) in the sum for the priority function 
#' G_W^D(s_i) used in the Fast Kulldorff method.
#' @inheritParams conditional_score_fun_EBP
conditional_score_fun_EBE <- function(c, b, q) {
  stop("Not yet implemented")
}