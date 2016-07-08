# Main function ----------------------------------------------------------------

#' Calculate the spatiotemporal ZIP scan statistic.
#' 
#' This function calculates the spatiotemporal zero-inflated Poisson scan 
#' statistic, along with the most likely cluster (MLC) and its Monte Carlo 
#' \eqn{p}-value. See details.
#' @param table A \code{data.table} with columns \code{duration, location, 
#'    count} and columns for distribution parameters that match the chosen type
#'    of scan statistic. See details.
#' @param zones A \code{set} of \code{set}s, the inner sets containing integer
#'    values that match those in the column \code{location} of the \code{table}
#'    argument of this function. See details.
#' @param n_replicates The number of Monte Carlo replications to perform in 
#'    order to calculate a p-value.
#' @param ... Arguments passed to internal functions. Arguments that can be
#'    passed here are \code{d}, the initial value for the excess zero indicators
#'    (default is 0.5), and \code{tol}, the threshold for the absolute 
#'    convergence criterion (default is 0.01).
#' @export
zip_scanstatistic <- function(table, zones, n_replicates, ...) {
  # input validation
  
  # Calculate statistics for observed data
  # extract max value
  maxdur <- table[, max(duration)]
  observed_statistics <- zip_calculations(table, zones, maxdur = maxdur, ...)
  scan_obs <- extract_scanstatistic(observed_statistics)
  
  replicate_scanstats <- zip_mcsim(table, 
                                   zones, 
                                   n_replicates, 
                                   maxdur = maxdur,
                                   ...)
  pval <- mc_pvalue(scan_obs, replicate_scanstats)
  
  list(data = table,
       zones = zones,
       n_replicates = n_replicates,
       replicates = replicate_scanstats,
       observed = observed_statistics,
       mlc = extract_mlc(observed_statistics),
       pvalue = pval)
}


# Simulation and hypothesis testing functions ----------------------------------

#' Randomly generate and add ZIP-distributed counts to a table.
#' 
#' This function randomly generates counts from a zero-inflated Poisson 
#' distribution according to the parameters on each row of the input 
#' \code{data.table}, and adds the counts to a new column \code{count}. 
#' @param table A \code{data.table} with columns \code{mean} and \code{p}. These
#'    correspond to the parameters \code{mu} and \code{sigma} in 
#'    \code{\link[gamlss.dist]{rZIP}}; the former is the Poisson mean and the 
#'    latter is the excess zero probability.
#' @param abs_tol A scalar; if the excess zero probability is less than this 
#'    value the count will be generated from a Poisson distribution.
#' @return The same table, with a new column \code{count}.
#' @importFrom gamlss.dist rZIP
#' @importFrom stats rpois
#' @keywords internal
generate_zip_counts <- function(table, abs_tol = 1e08) {
  # Note: rZIP returns a numeric vector
  table[, count := 0]
  table[p < abs_tol, count := as.numeric(rpois(.N, mean))]
  table[p >= abs_tol, count := gamlss.dist::rZIP(.N, mu = mean, sigma = p)][]
}

#' Simulate a single expectation-based ZIP-EM scan statistic.
#' 
#' Simulate zero-inflated -distributed data according to the supplied parameters 
#' and calculate the value of the expectation-based Poisson scan statistic.
#' @param table A \code{data.table} with columns \code{location, duration, p,
#'    mean}. The column \code{mean} contains the Poisson means, the column 
#'    \code{p} contains the excess zero probabilities.
#' @inheritParams partition_zones
#' @param ... Arguments passed to \code{\link{zip_calculations}}.
#' @return A scalar; the expectation-based ZIP-EM scan statistic for the 
#'    simulated data.
#' @importFrom magrittr %>%
#' @keywords internal
simulate_zip_scanstatistic <- function(table, zones, ...) {
  table[, .(p, mean), by = .(location, duration)] %>%
    generate_zip_counts %>% 
    zip_calculations(zones = zones, ...) %>%
    extract_scanstatistic
}

#' Monte Carlo simulation of expectation-based Poisson scan statistics.
#' 
#' This function generates \code{n_replicates} Poisson-distributed data sets 
#' according to the parameters in the input table, and calculates the value of
#' the scan statistic for each generated data set using the supplied 
#' \code{zones}.
#' @param table A \code{data.table} with columns \code{location, duration, p,
#'    mean}.
#' @inheritParams partition_zones
#' @param n_replicates A positive integer; the number of replicate scan 
#'    statistics to generate.
#' @param ... Arguments passed to \code{\link{simulate_zip_scanstatistic}}.
#' @return A numeric vector of length \code{n_replicates}.
#' @importFrom magrittr %>%
#' @keywords internal
zip_mcsim <- function(table, zones, n_replicates = 999L, ...) {
  replicate(n_replicates, simulate_zip_scanstatistic(table, zones, ...))
}


# Functions with data.table input ----------------------------------------------

#' Estimate the ZIP excess zero indicators under null hypothesis of no outbreak.
#' 
#' Estimate the excess zero indicator for a zero-inflated Poisson distribution,
#' under the null hypothesis of no outbreak. The estimate is added as a new 
#' column to the input table.
#' @param table A \code{data.table} with columns \code{p, mean, count}. \code{p}
#'    is the given/estimated probability of an excess zero, and \code{mean} is
#'    the estimated Poisson mean. \code{count} is the observed count.
#' @return The same table, \strong{modified} with a new column \code{ddagger}.
#' @keywords internal
estimate_d_dagger <- function(table) {
  table[, ddagger := estimate_d(p, mean, count)][]
}

#' Calculate the ZIP window statistic over all durations, for a given zone.
#' 
#' This function calculates the zero-inflated Poisson statistic for a given 
#' spatial zone, for all durations considered.
#' @param table A \code{data.table} with columns \code{duration, location, 
#'    p, mean, count}.
#' @param maxdur An integer; the maximum duration considered.
#' @param ... Arguments passed to \code{\link{window_zip_statistic}}.
#' @return A list with two elements:
#' \describe{
#'   \item{duration}{Vector of integers from 1 to \code{maxdur}.}
#'   \item{statistic}{Numeric vector containing the ZIP statistics corresponding
#'   to each duration, for the given spatial zone.}
#' }
#' @keywords internal
calc_zipstat_over_duration <- function(table, maxdur, ...) {
  stat <- rep(0, maxdur)
  for (t in seq(maxdur)) {
    stat[t] <- table[duration <= t, window_zip_statistic(p, mean, count, ...)]
  }
  list(duration = seq(maxdur), statistic = stat)
}

#' Calculates the ZIP statistic for each space-time window.
#' 
#' Calculates the zero-inflated Poisson statistic for each space-time window,
#' using the EM algorithm.
#' @param table A \code{data.table} with columns \code{zone, location, 
#'    duration, p, mean, count}.
#' @param ... Any of the following named parameters:
#' \describe{
#'   \item{maxdur}{As in \code{\link{calc_zipstat_over_duration}}.}
#'   \item{d_init}{As in \code{link{window_zip_statistic}}.}
#'   \item{tol}{As in \code{link{window_zip_statistic}}.}
#' }
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#' @keywords internal
zip_statistic <- function(table, ...) {
  table[, calc_zipstat_over_duration(.SD, ...), by = .(zone)]
}

#' Calculate the (logarithm of the) ZIP statistic for each space-time window.
#' 
#' Calculate the logarithm of the ZIP statistic for each space-time window, 
#' summing over the locations and times in the window.
#' @param table A \code{data.table} with columns \code{location, duration, mean,
#'    p, count}. The column \code{mean} contains the Poisson means, the column
#'    \code{p} contains the excess zero probabilities.
#' @inheritParams partition_zones
#' @param ... Arguments passed to \code{\link{simulate_zip_scanstatistic}}.
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#' @importFrom magrittr %>%
#' @keywords internal
zip_calculations <- function(table, zones, ...) {
  table %>%
    zone_joiner(zones = zones, keys = c("zone", "duration")) %>%
    zip_statistic(...)
}


# Functions with scalar/vector input -------------------------------------------

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
#' @keywords internal
estimate_zip_relrisk <- function(d, mu, y) {
  max(1, sum(y * (1 - d)) / sum(mu * (1 - d)))
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
#' @keywords internal
estimate_d <- function(p, mu, y) {
  d <- rep(0, length(y))
  zero <- y == 0L
  d[zero] <- p[zero] / (p[zero] + (1 - p[zero]) * exp(-mu[zero]))
  d
}

#' Calculate a term in the sum of the logarithm of the ZIP window statistic.
#' 
#' This function calculates a term which appears in the sum of the logarithm of
#' the zero-inflated Poisson statistic for a given space-time window.
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
#' @keywords internal
zip_statistic_term <- function(q, p, dstar, ddagger, mu, y) {
  zip_statistic_factor(p, dstar, q * mu, y) - 
    zip_statistic_factor(p, ddagger, mu, y)
}

#' Factor(s) in the product of the EB-ZIP window statistic.
#' 
#' Computes one or more factors in the product of the numerator or denominator
#' of the EB-ZIP window statistic (i.e. the one calculated for a given 
#' space-time window W). Should be able to handle cases where p = 0.
#' @param p Numeric vector of excess zero probabilities.
#' @param d Numeric vector of estimates of the excess zero indicators. Of same 
#'    length as \code{p}.
#' @param mu Numeric vector of given/estimated Poisson means. Of same length as 
#'    \code{p}.
#' @param y Integer vector of observed counts. Of same length as \code{p}.
#' @param tol Scalar; probability p below this is considered equal to zero.
#' @return A numeric vector of same length as input vector \code{p}.
#' @importFrom stats dpois
#' @keywords internal
zip_statistic_factor <- function(p, d, mu, y, tol = 1e-08) {
  res <- rep(0, length(p))
  p_is_zero <- p < tol
  res[p_is_zero] <- dpois(y[p_is_zero], mu[p_is_zero], log = TRUE)
  res[!p_is_zero] <- d[!p_is_zero] * log(p[!p_is_zero]) + 
    (1 - d[!p_is_zero]) * (log(1 - p[!p_is_zero]) + 
                             dpois(y[!p_is_zero], mu[!p_is_zero], log = TRUE))
  res
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
#' @keywords internal
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
#' @param ... Named parameters passed to \code{\link{zip_em_estimates}}.
#' @return A scalar, the (logarithm of the) ZIP statistic.
#' @keywords internal
window_zip_statistic <- function(p, mu, y, ...) {
  em <- zip_em_estimates(p, mu, y, ...)
  sum(zip_statistic_term(em$q, p, em$dstar, estimate_d(p, mu, y), mu, y))
}