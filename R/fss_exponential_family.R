
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
#'    and \code{b} (baseline) and returning a scalar.
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

# Gaussian ---------------------------------------------------------------------

#' Calculates the largest value for which the Gaussian score function is zero.
#' @param c A scalar; the observed value.
#' @param b A scalar; the expected value.
#' @return A scalar greater than or equal to 1.
#' @keywords internal
gaussian_qmax <- function(c, b) {
  ifelse(c > b, 2 * c / b - 1, 1)
}

#' Compute Gaussian priority function values cumulatively over time.
#' 
#' Given a matrices with observed and expected values (baselines) for each 
#' timepoint (row) and location/data stream (column), sum counts and baselines 
#' cumulatively backwards in time, and compute the Gaussian priority function
#' value for each column and subset of time.
#' @param counts A matrix of observed counts. Rows represent timepoints, ordered
#'    from most recent to most distant. Columns represent e.g. locations or
#'    data streams.
#' @param baselines A matrix of expected counts with the same dimensions as
#'    \code{counts}.
#' @param scalar_priority_fun A function taking two arguments \code{c} (count)
#'    and \code{b} (baseline) and returning a scalar.
#' @return A matrix with the same dimensions as \code{counts}. The \eqn{i}th 
#'    element in each column contains the Gaussian priority function value for 
#'    the window of time (duration) stretching from 1 to \eqn{i}, for that
#'    column.
#' @importFrom purrr map2_dbl
#' @keywords internal
gaussian_priority <- function(counts, 
                              baselines, 
                              scalar_priority_fun = gaussian_qmax) {
  agg_c <- as.vector(apply(counts, 2, cumsum))
  agg_b <- as.vector(apply(baselines, 2, cumsum))
  prios <- map2_dbl(agg_c, agg_b, scalar_priority_fun)
  matrix(prios, nrow(counts), ncol(counts))
}

#' Compute the score for an individual Gaussian observation.
#' @param c A scalar; the observed value.
#' @param b A scalar; the expected value.
#' @param s2 A scalar; the variance.
#' @return A scalar; the score.
#' @keywords internal
gaussian_lambda <- function(c, b, s2) {
  q <- c / b
  c * b * (q - 1) / s2 + b^2 * (1 - q^2) / (2 * s2)
}

#' Compute the Gaussian score for each priority subset, for each duration.
#' @inheritParams gaussian_priority
#' @param variances A matrix of variances with the same dimensions as 
#'    \code{counts}.
#' @param priority_indices A matrix of the same size as the input. On each row
#'    (duration), column numbers are given in order of priority.
#' @return A matrix of the same dimension as the input matrices.
#' @importFrom purrr pmap_dbl
#' @keywords internal
gaussian_score <- function(counts, baselines, variances, priority_indices) {
  # Since the sum of Gaussian random variables is Gaussian with parameters 
  # (expected value and variance) equal to the sum of the individual parameters, 
  # compute cumulative sums over duration (rows) for each column. Then 
  # re-arrange values in each row by priority, and use the same property again 
  # by computing cumulative sums over locations/data streams (columns) for each 
  # row.
  counts <- sum_reorder_sum(counts, priority_indices)
  baselines <- sum_reorder_sum(baselines, priority_indices)
  variances <- sum_reorder_sum(variances, priority_indices)
  args <- list(c = counts, b = baselines, s2 = variances)
  
  # Compute scores for corresponding element pairs in the count and baseline
  # matrices
  matrix(pmap_dbl(.l = args, .f = gaussian_lambda), 
         nrow(priority_indices), 
         ncol(priority_indices))
}

# Exponential ---------------------------------------------------------------------

#' Calculates the largest value for which the exponential score function is
#' zero.
#' @param c A scalar; the observed value.
#' @param b A scalar; the expected value.
#' @return A scalar greater than or equal to 1.
#' @keywords internal
exponential_qmax <- function(c, b) {
  -c/b * 1 / emdbook::lambertW(-c/b * exp(-c/b), b = 0)
}

#' Compute exponential priority function values cumulatively over time.
#'
#' Given a matrices with observed and expected values (baselines) for each
#' timepoint (row) and location/data stream (column), sum counts and baselines
#' cumulatively backwards in time, and compute the exponential priority function
#' value for each column and subset of time.
#' @param counts A matrix of observed counts. Rows represent timepoints, ordered
#'    from most recent to most distant. Columns represent e.g. locations or
#'    data streams.
#' @param baselines A matrix of expected counts with the same dimensions as
#'    \code{counts}.
#' @param scalar_priority_fun A function taking two arguments \code{c} (count)
#'    and \code{b} (baseline) and returning a scalar.
#' @return A matrix with the same dimensions as \code{counts}. The \eqn{i}th
#'    element in each column contains the exponential priority function value
#'    for the window of time (duration) stretching from 1 to \eqn{i}, for that
#'    column.
#' @importFrom purrr map2_dbl
#' @keywords internal
exponential_priority <- function(counts,
                                 baselines,
                                 scalar_priority_fun = exponential_qmax) {
  agg_c <- as.vector(apply(counts, 2, cumsum))
  agg_b <- as.vector(apply(baselines, 2, cumsum))
  prios <- map2_dbl(agg_c, agg_b, scalar_priority_fun)
  matrix(prios, nrow(counts), ncol(counts))
}

#' Compute the score for an individual exponential observation.
#' @param c A scalar; the observed value.
#' @param b A scalar; the expected value.
#' @return A scalar; the score.
#' @keywords internal
exponential_lambda <- function(c, b) {
  q <- c / b
  q * (1 - 1 / q) - log(q)
}

#' Compute the exponential score for each priority subset, for each duration.
#' @inheritParams exponential_priority
#' @param priority_indices A matrix of the same size as the input. On each row
#'    (duration), column numbers are given in order of priority.
#' @return A matrix of the same dimension as the input matrices.
#' @importFrom purrr pmap_dbl
#' @keywords internal
exponential_score <- function(counts, baselines, priority_indices) {

  counts <- sum_reorder_sum(counts, priority_indices)
  baselines <- sum_reorder_sum(baselines, priority_indices)
  variances <- sum_reorder_sum(variances, priority_indices)
  args <- list(c = counts, b = baselines)

  # Compute scores for corresponding element pairs in the count and baseline
  # matrices
  matrix(pmap_dbl(.l = args, .f = exponential_lambda),
         nrow(priority_indices),
         ncol(priority_indices))
}

