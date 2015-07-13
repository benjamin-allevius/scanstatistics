

#' Calculate the (logarithm of the) ZIP statistic for each space-time window.
#' 
#' Calculate the logarithm of the ZIP statistic for each space-time window, 
#' summing over the locations and times in the window.
#' @param table A \code{data.table} 
# calculate_zip_statistic <- function(table, region) {
#   table %>%
#     region_joiner(regions = regions, keys = c("region", "duration")) %>%
#     # calculateq for each window, then the location-time (individual obs) parts of statistic sum, then statistic value 
#     region_sum(sumcols = c("count", "mean")) %>%
#     cumsum_duration(sumcols = c("count", "mean"), bycols = c("region")) %>%
#     poisson_relrisk %>%
#     poisson_statistic
# #   table[, 
# #     statistic := sum(zip_stat_formula(relrisk, p, dstar, ddagger, mean, count)),
# #     by = .(region, duration)]
# }


#' Formula for a term in the sum of the logarithm of the ZIP window statistic.
#' @param q Scalar; the relative risk.
#' @param p Numeric vector of excess zero probabilities.
#' @param dstar Numeric vector of estimates of the excess zero indicators, under 
#'    the alternative hypothesis of an outbreak. Of same length as \code{p}.
#' @param ddagger Numeric vector of estimates of the excess zero indicators, 
#'    under the null hypothesis of no outbreak. Of same length as \code{p}.
#' @param mu Numeric vector of given/estimated Poisson means. Of same length as 
#'    \code{p}.
#' @param y Integer vector of observed counts. Of same length as \code{p}.
#' @return A numeric vector of same length as input vector \code{p}.
zip_statistic_term <- function(q, p, dstar, ddagger, mu, y) {
  (dstar - ddagger) * (log(p) - log(1 - p) - y * log(mu) + lfactorial(y)) +
    (1 - dstar) * (y * log(q) - q * mu) + (1 - ddagger) * mu
}


#' Estimate the indicators of excess zeros for a ZIP distribution.
#' 
#' Given counts and (estimated) Poisson means and excess zero probabilities, 
#' this function estimates the indicator variable for an excess zero, for each
#' count.
#' @param p A numeric vector, each element being the given/estimated 
#'    probability of an excess zero for the corresponding count in \code{y}.
#' @param mu A numeric vector, each element being the given/estimated Poisson 
#'    meanfor the corresponding count in \code{y}. Of same length as \code{p}.
#' @param y An integer vector containing the observed counts. Of same length as 
#'    \code{p}.
#' @return A numeric vector, of same length as the input vector \code{p}.
estimate_d <- function(p, mu, y) {
  res <- rep(0, length(y))
  zero <- y == 0L
  res[zero] <- p[zero] / (p[zero] + (1 - p[zero]) * exp(-mu[zero]))
  res
}


#' Estimate the ZIP excess zero indicators under null hypothesis of no outbreak.
#' 
#' Estimate the excess zero indicator for a zero-inflated Poisson distribution,
#' under the null hypothesis of no outbreak. The estimate is added as a new 
#' column to the input table.
#' @param table A \code{data.table} with columns \code{p, mean, count}. \code{p}
#'    is the given/estimated probability of an excess zero, and \code{mean} is
#'    the estimated Poisson mean. \code{count} is the observed count.
#' @return The same table, \strong{modified} with a new column \code{ddagger}.
estimate_d_dagger <- function(table) {
  table[, ddagger := estimate_d(p, mean, count)][]
}

#' Estimate the relative risk for an outbreak using a ZIP distribution.
#' 
#' Scalar estimate of the relative risk \eqn{q} for zero-inflated Poisson data.
#' @param d A vector of indicator variables for whether the corresponding count
#'    in the argument \code{y} is an excess zero or not. Can also be estimates 
#'    of these indicators.
#' @param mu A vector of given/estimated Poisson means, of same length as 
#'    \code{d}.
#' @param y An integer vector of the observed counts, of same length as 
#'    \code{d}.
#' @return A scalar, the estimated relative risk.
estimate_zip_relrisk <- function(d, mu, y) {
  max(1, sum(y * (1 - d)) / sum(mu * (1 - d)))
}

#' Estimates the ZIP relative risk and excess zero indicators for a window.
#' 
#' For a single spatial or space-time window, this function uses the EM 
#' algorithm to estimate the relative risk and the excess zero indicators for 
#' counts assumed to be generated from a zero-inflated Poisson distribution.
#' @param p A numeric vector of the given/estimated excess zero probabilities 
#'    corresponding to each count.
#' @param mu A numeric vector of the given/estimated Poisson means corresponding 
#'    to each count. Of same length as \code{p}.
#' @param y An integer vector of the observed counts, of same length as 
#'    \code{p}.
#' @param d_init A scalar between 0 and 1. The initial guess for the estimate of 
#'    the excess zero indicator.
#' @param tol A scalar between 0 and 1. It is the absolute tolerance criterion
#'    for the estimate of the excess zero indicator; convergence is reached when
#'    two successive elements in the sequence of estimates have an absolute 
#'    difference less than \code{tol}.
#' @return A list with two elements:
#' \describe{
#'   \item{q}{Scalar estimate of the relative risk.}
#'   \item{dstar}{Estimates of the excess zero indicator variables.}
#' }
zip_em_estimates <- function(p, mu, y, d_init = 0.5, tol = 0.01) {
  d_prev <- ifelse(y > 0, rep(0, length(y)), rep(d_init, length(y)))
  q <- estimate_zip_relrisk(d_prev, mu, y)
  d <- estimate_d(p, q * mu, y)
  
  left <- abs(d - d_prev) >= tol
  while (any(left)) {
    d_prev[left] <- d[left]
    q <- estimate_zip_relrisk(d, mu, y)
    d[left] <- estimate_d(p[left], q * mu[left], y[left])
    left[left] <- abs(d[left] - d_prev[left]) >= tol
  }
  list(q = estimate_zip_relrisk(d, mu, y), dstar = d)
}

#' Calculate the ZIP statistic for a single space-time window.
#' 
#' Calculate the single-window statistic for the zero-inflated Poisson 
#' distribution using the EM algorithm.
#' @inheritParams zip_em_estimates
#' @return A scalar, the (logarithm of the) ZIP statistic.
window_zip_statistic <- function(p, mu, y, d_init = 0.5, tol = 0.01) {
  em <- zip_em_estimates(p, mu, y, d_init = 0.5, tol = 0.01)
  sum(zip_statistic_term(em$q, p, em$dstar, estimate_d(p, mu, y), mu, y))
}

#' Calculate the ZIP window statistic over all durations, for a given region.
#' @param table A \code{data.table} with columns \code{duration, location, 
#'    p, mean, count}.
#' @param maxdur An integer; the maximum duration considered.
#' @return A list with two elements:
#' \describe{
#'   \item{duration}{Vector of integers from 1 to \code{maxdur}.}
#'   \item{statistic}{Numeric vector containing the ZIP statistics corresponding
#'   to each duration, for the given spatial zone.}
#' }
calc_zipstat_over_duration <- function(table, maxdur) {
  stat <- rep(0, maxdur)
  for (t in seq(maxdur)) {
    stat[t] <- table[duration <= t, window_zip_statistic(p, mean, count)]
  }
  list(duration = seq(maxdur), statistic = stat)
}

zip_statistic <- function(table, ...) {
  table[, calc_zipstat_over_duration(.SD, ...), by = .(region)]
}