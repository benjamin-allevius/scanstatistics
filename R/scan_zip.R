# Functions in this file:
#   scan_zip
#   generate_zip_counts
#   simulate_zip_scanstatistic
#   zip_mcsim
#   estimate_d_dagger
#   calc_zipstat_over_duration
#   zip_statistic
#   zip_calculations
#   estimate_zip_relrisk
#   estimate_d
#   zip_statistic_term
#   zip_statistic_factor
#   zip_em_estimates
#   window_zip_statistic

# Main function ----------------------------------------------------------------

#' Calculate the ZIP scan statistic.
#' 
#' Calculate the expectation-based zero-inflated Poisson scan statistic by 
#' supplying a \code{data.table} of observed counts and pre-computed expected 
#' values and structural zero probabilities for each location and time. A 
#' p-value for the observed scan statistic can be obtained by Monte Carlo 
#' simulation.
#' 
#' @param table A \code{data.table} with columns 
#'    \code{location, duration, count, mu, p}. The \code{location} column 
#'    should consist of integers that are unique to each location. The 
#'    \code{duration} column should also consist of integers, starting at 1 for 
#'    the most recent time period and increasing in reverse chronological order.
#'    The column \code{mu} should contain the estimated Poisson expected value 
#'    parameters, and the column \code{p} the estimated structural zero 
#'    probabilities.
#' @param zones A \code{set} of zones, each zone itself a set containing one or 
#'    more locations of those found in \code{table}.
#' @param n_mcsim A non-negative integer; the number of replicate scan 
#'    statistics to generate in order to calculate a p-value.
#' @param ... Arguments passed to internal functions. Arguments that can be
#'    passed here are \code{d}, the initial value for the excess zero indicators
#'    (default is 0.5), and \code{tol}, the threshold for the absolute 
#'    convergence criterion (default is 0.01).
#' @return An object of class \code{scanstatistics}. It has the following 
#'    fields:
#'    \describe{
#'     \item{observed}{A \code{data.table} containing the value of the 
#'                     statistic calculated for each zone-duration combination,
#'                     for the observed data. The scan statistic is the maximum
#'                     value of these calculated statistics.}
#'     \item{replicated}{A numeric vector of length \code{n_mcsim} containing 
#'                       the values of the scanstatistics calculated by Monte
#'                       Carlo simulation.}
#'     \item{mlc}{A \code{data.table} containing the zone, duration, and 
#'                scanstatistic.}
#'     \item{pvalue}{The p-value calculated from Monte Carlo replications.}
#'     \item{distribution}{The assumed distribution of the data; "zero-inflated 
#'                         Poisson" in this case.}
#'     \item{type}{The type of scan statistic; "Expectation-based" in this 
#'                 case.}
#'     \item{zones}{The set of zones that was passed to the function as input.}
#'     \item{n_locations}{The number of locations in the data.}
#'     \item{n_zones}{The number of zones.}
#'     \item{max_duration}{The maximum anomaly duration considered.}
#'    }
#' @export
#' @concept zero-inflated poisson scanstatistic
#' @details For the expectation-based zero-inflated Poisson scan statistic 
#'    (Kjellson 2015), the null hypothesis of no anomaly holds that the count 
#'    observed at each location \eqn{i} and duration \eqn{t} (the number of time 
#'    periods before present) has a zero-inflated Poisson distribution with 
#'    expected value parameter \eqn{\mu_{it}} and structural zero probability
#'    \eqn{p_{it}}:
#'    \deqn{
#'      H_0 : Y_{it} \sim \textrm{ZIP}(\mu_{it}, p_{it}).
#'    }
#'    This holds for all locations \eqn{i = 1, \ldots, m} and all durations 
#'    \eqn{t = 1, \ldots,T}, with \eqn{T} being the maximum duration considered.
#'    Under the alternative hypothesis, there is a space-time window \eqn{W}
#'    consisting of a spatial zone \eqn{Z \subset \{1, \ldots, m\}} and a time 
#'    window \eqn{D \subseteq \{1, \ldots, T\}} such that the counts in that
#'    window have their Poisson expected value parameters inflated by a factor 
#'    \eqn{q_W > 1} compared to the null hypothesis:
#'    \deqn{
#'    H_1 : Y_{it} \sim \textrm{ZIP}(q_W \mu_{it}, p_{it}), ~~(i,t) \in W.
#'    }
#'    For locations and durations outside of this window, counts are assumed to
#'    be distributed as under the null hypothesis. The sets \eqn{Z} considered 
#'    are those specified in the argument \code{zones}, while the maximum 
#'    duration \eqn{T} is taken as the maximum value in the column 
#'    \code{duration} of the input \code{table}. 
#'    
#'    For each space-time window \eqn{W} considered, (the log of) a likelihood 
#'    ratio is computed using the distributions under the alternative and null 
#'    hypotheses, and the expectation-based Poisson scan statistic is calculated 
#'    as the maximum of these quantities over all space-time windows. The 
#'    expectation-maximization (EM) algorithm is used to obtain maximum 
#'    likelihood estimates. Point estimates of the parameters \eqn{\mu_{it}} 
#'    must be specified in the column \code{mu} of the argument \code{table} 
#'    before this function is called.
#' @references 
#'    Kjellson, B. (2015), \emph{Spatiotemporal Outbreak Detection: A Scan 
#'    Statistic Based on the Zero-Inflated Poisson Distribution}, (Master 
#'    Thesis, Stockholm University),
#'    \href{http://goo.gl/6Q89ML}{Link to PDF}.
#' @examples 
#' # Simple example
#' set.seed(1)
#' table <- scanstatistics:::create_table(list(location = 1:4, duration = 1:4),
#'                                         keys = c("location", "duration"))
#' table[, mu := 3 * location]
#' table[, p := runif(.N, 0, 0.3)]
#' table[, count := gamlss.dist::rZIP(.N, mu = mu, sigma = p)]
#' table[location %in% c(1, 4) & duration < 3, 
#'       count := gamlss.dist::rZIP(.N, mu = 2 * mu, sigma = p)]
#' zones <- scanstatistics:::powerset_zones(4)
#' result <- scan_poisson(table, zones, 100)
#' result
scan_zip <- function(table, zones, n_mcsim = 0, ...) {
  validate_scan(table, zones, c("count", "mu", "duration", "location", "p"))
  maxdur <- table[, max(duration)]
  scanstatistic_object(zip_calculations(table, zones, maxdur, ...), 
                       zip_mcsim(table, zones, n_mcsim, maxdur, ...),
                       list(table = table,
                            zones = zones, 
                            distribution = "zero-inflated Poisson",
                            type = "Expectation-based"))
}


# Simulation and hypothesis testing functions ----------------------------------

#' Randomly generate and add ZIP-distributed counts to a table.
#' 
#' This function randomly generates counts from a zero-inflated Poisson 
#' distribution according to the parameters on each row of the input 
#' \code{data.table}, and adds the counts to a new column \code{count}. 
#' @param table A \code{data.table} with columns \code{mu} and \code{p}. These
#'    correspond to the parameters \code{mu} and \code{sigma} in 
#'    \code{\link[gamlss.dist]{rZIP}}; the former is the Poisson expected value
#'    parameter and the latter is the excess zero probability.
#' @param abs_tol A scalar; if the excess zero probability is less than this 
#'    value the count will be generated from a Poisson distribution.
#' @return The same table, with a new column \code{count}.
#' @importFrom gamlss.dist rZIP
#' @importFrom stats rpois
#' @keywords internal
generate_zip_counts <- function(table, abs_tol = 1e-03) {
  # Note: rZIP returns a numeric vector
  table[, count := 0]
  table[p < abs_tol, count := as.numeric(rpois(.N, mu))]
  table[p >= abs_tol, count := gamlss.dist::rZIP(.N, mu = mu, sigma = p)][]
}

#' Simulate a single expectation-based ZIP-EM scan statistic.
#' 
#' Simulate zero-inflated -distributed data according to the supplied parameters 
#' and calculate the value of the expectation-based Poisson scan statistic.
#' @param table A \code{data.table} with columns \code{location, duration, p,
#'    mu}. The column \code{mu} contains the Poisson expected values, the column 
#'    \code{p} contains the excess zero probabilities.
#' @inheritParams partition_zones
#' @param ... Arguments passed to \code{\link{zip_calculations}}.
#' @return A scalar; the expectation-based ZIP-EM scan statistic for the 
#'    simulated data.
#' @importFrom magrittr %>%
#' @keywords internal
simulate_zip_scanstatistic <- function(table, zones, ...) {
  table[, .(p, mu), by = .(location, duration)] %>%
    generate_zip_counts %>% 
    zip_calculations(zones = zones, ...) %>%
    extract_scanstatistic
}

#' Monte Carlo simulation of expectation-based Poisson scan statistics.
#' 
#' This function generates \code{n_mcsim} Poisson-distributed data sets 
#' according to the parameters in the input table, and calculates the value of
#' the scan statistic for each generated data set using the supplied 
#' \code{zones}.
#' @param table A \code{data.table} with columns \code{location, duration, p,
#'    mu}.
#' @inheritParams partition_zones
#' @param n_mcsim A positive integer; the number of replicate scan 
#'    statistics to generate.
#' @param ... Arguments passed to \code{\link{simulate_zip_scanstatistic}}.
#' @return A numeric vector of length \code{n_mcsim}.
#' @importFrom magrittr %>%
#' @keywords internal
zip_mcsim <- function(table, zones, n_mcsim = 0, ...) {
  if (n_mcsim > 0) {
    return(replicate(n_mcsim, simulate_zip_scanstatistic(table, zones, ...)))
  } else {
    return(numeric(0))
  }
}


# Functions with data.table input ----------------------------------------------

#' Estimate the ZIP excess zero indicators under null hypothesis of no outbreak.
#' 
#' Estimate the excess zero indicator for a zero-inflated Poisson distribution,
#' under the null hypothesis of no outbreak. The estimate is added as a new 
#' column to the input table.
#' @param table A \code{data.table} with columns \code{p, mu, count}. \code{p}
#'    is the given/estimated probability of an excess zero, and \code{mu} is
#'    the estimated Poisson expected value parameter. \code{count} is the 
#'    observed count.
#' @return The same table, \strong{modified} with a new column \code{ddagger}.
#' @keywords internal
estimate_d_dagger <- function(table) {
  table[, ddagger := estimate_d(p, mu, count)][]
}

#' Calculate the ZIP window statistic over all durations, for a given zone.
#' 
#' This function calculates the zero-inflated Poisson statistic for a given 
#' spatial zone, for all durations considered.
#' @param table A \code{data.table} with columns \code{duration, location, p, 
#'    mu, count}.
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
    stat[t] <- table[duration <= t, window_zip_statistic(p, mu, count, ...)]
  }
  list(duration = seq(maxdur), statistic = stat)
}

#' Calculates the ZIP statistic for each space-time window.
#' 
#' Calculates the zero-inflated Poisson statistic for each space-time window,
#' using the EM algorithm.
#' @param table A \code{data.table} with columns \code{zone, location, duration, 
#'    p, mu, count}.
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
#' @param table A \code{data.table} with columns \code{location, duration, mu,
#'    p, count}. The column \code{mu} contains the Poisson expected value 
#'    parameters, the column \code{p} contains the excess zero probabilities.
#' @inheritParams partition_zones
#' @param ... Arguments passed to \code{\link{simulate_zip_scanstatistic}}.
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#' @importFrom magrittr %>%
#' @keywords internal
zip_calculations <- function(table, zones, ...) {
  table[, .(location, duration, count, mu, p)] %>%
    join_zones(zones = zones, keys = c("zone", "duration")) %>%
    zip_statistic(...)
}


# Functions with scalar/vector input -------------------------------------------

#' Estimate the relative risk for an outbreak using a ZIP distribution.
#' 
#' Scalar estimate of the relative risk \eqn{q} for zero-inflated Poisson data.
#' @param d A vector of indicator variables for whether the corresponding count
#'    in the argument \code{y} is an excess zero or not. Can also be estimates 
#'    of these indicators.
#' @param mu A vector of given/estimated Poisson expected value parameters, of 
#'    same length as \code{d}.
#' @param y An integer vector of the observed counts, of same length as 
#'    \code{d}.
#' @return A scalar, the estimated relative risk.
#' @keywords internal
estimate_zip_relrisk <- function(d, mu, y) {
  max(1, sum(y * (1 - d)) / sum(mu * (1 - d)))
}

#' Estimate the indicators of excess zeros for a ZIP distribution.
#' 
#' Given counts and (estimated) Poisson expected value parameters and excess 
#' zero probabilities, this function estimates the indicator variable for an 
#' excess zero, for each count.
#' @param p A numeric vector, each element being the given/estimated 
#'    probability of an excess zero for the corresponding count in \code{y}.
#' @param mu A numeric vector, each element being the given/estimated Poisson 
#'    expected value parameter for the corresponding count in \code{y}. Of same 
#'    length as \code{p}.
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
#' @param mu Numeric vector of given/estimated Poisson expected value 
#'    parameters. Of same length as \code{p}.
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
#' @param mu Numeric vector of given/estimated Poisson expected value 
#' parameters. Of same length as \code{p}.
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
#' @param mu A numeric vector of the given/estimated Poisson expected value
#'    parameters corresponding to each count. Of same length as \code{p}.
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