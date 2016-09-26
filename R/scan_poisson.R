# Functions in this file:
#   scan_poisson
#   generate_poisson_counts
#   simulate_poisson_scanstatistic
#   poisson_mcsim
#   poisson_calculations
#   poisson_statistic
#   poisson_relrisk

# Main function ----------------------------------------------------------------

#' Calculate the Poisson scan statistic.
#' 
#' Calculate the expectation-based Poisson scan statistic by supplying a 
#' \code{data.table} of observed counts and pre-computed expected value 
#' parameters for each location and time. A p-value for the observed scan
#' statistic can be obtained by Monte Carlo simulation.
#' 
#' @param table A \code{data.table} with columns 
#'    \code{location, duration, count, mu}. The \code{location} column should 
#'    consist of integers that are unique to each location. The 
#'    \code{duration} column should also consist of integers, starting at 1 for 
#'    the most recent time period and increasing in reverse chronological order.
#'    The column \code{mu} should contain the estimated Poisson expected value
#'    parameter.
#' @param zones A \code{set} of zones, each zone itself a 
#'    set containing one or more locations of those found in \code{table}.
#' @param n_mcsim A non-negative integer; the number of replicate scan 
#'    statistics to generate in order to calculate a p-value.
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
#'     \item{distribution}{The assumed distribution of the data; "Poisson" in
#'                         this case.}
#'     \item{type}{The type of scan statistic; "Expectation-based" in this 
#'                 case.}
#'     \item{zones}{The set of zones that was passed to the function as input.}
#'     \item{n_locations}{The number of locations in the data.}
#'     \item{n_zones}{The number of zones.}
#'     \item{max_duration}{The maximum anomaly duration considered.}
#'    }
#' @export
#' @concept poisson scanstatistic
#' @details For the expectation-based Poisson scan statistic, the null 
#'    hypothesis of no anomaly holds that the count observed at each location 
#'    \eqn{i} and duration \eqn{t} (the number of time periods before present) 
#'    is Poisson-distributed with expected value \eqn{\mu_{it}}:
#'    \deqn{
#'      H_0 : Y_{it} \sim \textrm{Poisson}(\mu_{it}),
#'    }
#'    for all locations \eqn{i = 1, \ldots, m} and all durations \eqn{t = 1,
#'    \ldots,T}, with \eqn{T} being the maximum duration considered.
#'    Under the alternative hypothesis, there is a space-time window \eqn{W}
#'    consisting of a spatial zone \eqn{Z \subset \{1, \ldots, m\}} and a time 
#'    window \eqn{D \subseteq \{1, \ldots, T\}} such that the counts in that
#'    window have their expected values inflated by a factor \eqn{q_W > 1} 
#'    compared to the null hypothesis:
#'    \deqn{
#'    H_1 : Y_{it} \sim \textrm{Poisson}(q_W \mu_{it}), ~~(i,t) \in W.
#'    }
#'    For locations and durations outside of this window, counts are assumed to
#'    be distributed as under the null hypothesis. The sets \eqn{Z} considered 
#'    are those specified in the argument \code{zones}, while the maximum 
#'    duration \eqn{T} is taken as the maximum value in the column 
#'    \code{duration} of the input \code{table}. For each space-time window
#'    \eqn{W} considered, (the log of) a likelihood ratio is computed using the 
#'    distributions under the alternative and null hypotheses, and the 
#'    expectation-based Poisson scan statistic is calculated as the maximum of 
#'    these quantities over all space-time windows.
#'    Point estimates of the parameters \eqn{\mu_{it}} must be specified in the
#'    column \code{mu} of the argument \code{table} before this function is 
#'    called.
#' @examples
#' # Simple example
#' set.seed(1)
#' table <- scanstatistics:::create_table(list(location = 1:4, duration = 1:4), 
#'                                         keys = c("location", "duration"))
#' table[, mu := 3 * location]
#' table[, count := rpois(.N, mu)]
#' table[location %in% c(1, 4) & duration < 3, count := rpois(.N, 2 * mu)]
#' zones <- scanstatistics:::powerset_zones(4)
#' result <- scan_poisson(table, zones, 100)
#' result
scan_poisson <- function(table, zones, n_mcsim = 0) {
  validate_scan(table, zones, c("count", "mu", "duration", "location"))
  scanstatistic_object(poisson_calculations(table, zones), 
                       poisson_mcsim(table, zones, n_mcsim),
                       list(table = table,
                            zones = zones, 
                            distribution = "Poisson",
                            type = "Expectation-based"))
}

# Simulation and hypothesis testing functions ----------------------------------

#' Randomly generate and add Poisson counts to a table.
#' 
#' This function randomly generates counts from a Poisson distribution according 
#' to the parameters on each row of the input \code{data.table}, and adds the 
#' counts to a new column \code{count}. 
#' @param table A \code{data.table} with at least the column \code{mu}, 
#'    corresponding to the parameter \code{lambda} in 
#'    \code{\link[stats]{rpois}}.
#' @importFrom stats rpois
#' @keywords internal
generate_poisson_counts <- function(table) {
  table[, count := rpois(.N, lambda = mu)][]
}

#' Simulate a single expectation-based Poisson scan statistic.
#' 
#' Simulate Poisson-distributed data according to the supplied parameters and
#' calculate the value of the expectation-based Poisson scan statistic.
#' @param table A \code{data.table} with columns \code{location, duration, mu}.
#' @inheritParams partition_zones
#' @return A scalar; the scan statistic for the simulated data.
#' @importFrom magrittr %>%
#' @keywords internal
simulate_poisson_scanstatistic <- function(table, zones) {
  table[, .(mu), by = .(location, duration)] %>%
    generate_poisson_counts %>% 
    poisson_calculations(zones = zones) %>%
    extract_scanstatistic
}

#' Monte Carlo simulation of expectation-based Poisson scan statistics.
#' 
#' This function generates \code{n_mcsim} Poisson-distributed data sets 
#' according to the parameters in the input table, and calculates the value of
#' the scan statistic for each generated data set using the supplied 
#' \code{zones}.
#' @inheritParams simulate_poisson_scanstatistic
#' @inheritParams partition_zones
#' @param n_mcsim A positive integer; the number of replicate scan 
#'    statistics to generate.
#' @return A numeric vector of length \code{n_mcsim}.
#' @keywords internal
poisson_mcsim <- function(table, zones, n_mcsim = 999L) {
  if (n_mcsim > 0) {
    return(replicate(n_mcsim, simulate_poisson_scanstatistic(table, zones)))
  } else {
    return(numeric(0))
  }
}


# Functions with data.table input ----------------------------------------------

#' Calculate the expectation-based Poisson statistic for each space-time window.
#' 
#' Calculate the expectation-based Poisson statistic for each space-time window,
#' given the initial data of counts and mus.
#' @param table A \code{data.table} with columns \code{location, duration, 
#'    count, mu}.
#' @param zones A \code{list} or \code{set} of zones, each zone itself a 
#'    set containing one or more locations of those found in \code{table}.
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#'    The column \code{statistic} contains the logarithm of the statistic for 
#'    each zone-duration combination.
#' @importFrom magrittr %>%
#' @keywords internal
poisson_calculations <- function(table, zones) {
  table[, .(location, duration, count, mu)] %>% 
    join_zones(zones = zones, keys = c("zone", "duration")) %>%
    zone_sum(sumcols = c("count", "mu")) %>%
    cumsum_duration(sumcols = c("count", "mu"), bycols = c("zone")) %>%
    poisson_relrisk %>%
    poisson_statistic
}

#' Calculate the expectation-based Poisson statistic for each space-time window.
#' 
#' This function calculates the logarithm of the expectation-based Poisson 
#' statistic for each space-time window, given already calculated relative risks
#' and aggregate counts and mus.
#' @param table A \code{data.table} with columns \code{zone, duration, count,
#'    mu, relrisk}. The columns \code{count} and \code{mu} contain the sums
#'    of the counts and mus for the locations inside each zone and up to the
#'    given duration. The column \code{relrisk} contains the maximum likelihood
#'    estimate for the relative risk for each zone-duration combination.
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#'    The column \code{statistic} contains the logarithm of the statistic for 
#'    each zone-duration combination.
#' @keywords internal
poisson_statistic <- function(table) {
  table[, .(statistic = log(relrisk) * count - (relrisk - 1) * mu),
        by = .(zone, duration)]
}

#' Add column for relative risk MLE to table with aggregate counts and means.
#' 
#' This function adds a column for the maximum likelihood estimate for the
#' relative risk (assumed to be 1 or greater). The table should contain the
#' aggregates (sums) of the counts and expected value parameters for each zone 
#' and duration.
#' @param table A \code{data.table} with columns \code{zone, duration, count,
#'    mu}. The latter two are the sums of observed counts and estimated expected
#'    value parameters for the locations comprising the zone and up to the given 
#'    duration (i.e. cumulative sum for the duration).
#' @return The same table, with an extra column \code{relrisk}.
#' @keywords internal
poisson_relrisk <- function(table) {
  table[, relrisk := max(1, count / mu), by = .(zone, duration)]
}