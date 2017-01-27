
# Poisson ----------------------------------------------------------------------

#' Calculates the largest value for which the Poisson score function is zero.
#' @param c A non-negative integer; the observed count.
#' @param b A positive scalar; the expected count.
#' @return A scalar greater than or equal to 1.
#' @importFrom emdbook lambertW
#' @keywords internal
poisson_qmax <- function(c, b) {
  r <- c / b
  ifelse(r > 1,
         -r * lambertW(-(1/r) * exp(-1/r), b = -1),
         1)
}

#' Compute Poisson priority function values cumulatively over time.
#' 
#' Given a matrices with observed and expected counts (baselines) for each 
#' timepoint (row) and location/data stream (column), sum counts and baselines 
#' cumulatively backwards in time, and compute the Poisson priority function
#' value for each column and subset of time.
#' @param counts A matrix of observed counts. Rows represent timepoints, ordered
#'    from most recent to most distant. Columns represent e.g. locations or
#'    data streams.
#' @param baselines A matrix of expected counts with the same dimensions as
#'    \code{counts}.
#' @param scalar_priority_fun A function taking two arguments \code{c} (count)
#'    and \code{b} (baseline) and 
#' @return A matrix with the same dimensions as \code{counts}. The \eqn{i}th 
#'    element in each column contains the Poisson priority function value for 
#'    the window of time (duration) stretching from 1 to \eqn{i}, for that
#'    column.
#' @importFrom purrr map2_dbl
#' @keywords internal
poisson_priority <- function(counts, 
                             baselines, 
                             scalar_priority_fun = poisson_qmax) {
  agg_c <- as.vector(apply(counts, 2, cumsum))
  agg_b <- as.vector(apply(baselines, 2, cumsum))
  prios <- suppressWarnings(map2_dbl(agg_c, agg_b, scalar_priority_fun))
  matrix(prios, nrow(counts), ncol(counts))
}

#' Compute the score for an individual Poisson observation.
#' @param c A non-negative integer; the observed count.
#' @param b A positive scalar; the expected count.
#' @return A scalar; the score.
#' @keywords internal
poisson_lambda <- function(c, b) {
  q <- c / b
  c * log(q) + b * (1 - q)
}

#' Compute the Poisson score for each priority subset, for each duration.
#' @inheritParams poisson_priority
#' @param priority_indices A matrix of the same size as the input. On each row
#'    (duration), column numbers are given in order of priority.
#' @return A matrix of the same dimension as the input matrices.
#' @importFrom purrr map2_dbl
#' @keywords internal
poisson_score <- function(counts, baselines, priority_indices) {
  # Since the sum of Poisson random variables is Poisson with parameter equal
  # to the sum of the individual parameters, compute cumulative sums over 
  # duration (rows) for each column. Then re-arrange values in each row by 
  # priority, and use the same property again by computing cumulative sums over
  # locations/data streams (columns) for each row.
  counts <- sum_reorder_sum(counts, priority_indices)
  baselines <- sum_reorder_sum(baselines, priority_indices)
  
  # Compute scores for corresponding element pairs in the count and baseline
  # matrices
  matrix(map2_dbl(counts, baselines, poisson_lambda), 
         nrow(priority_indices), 
         ncol(priority_indices))
}
