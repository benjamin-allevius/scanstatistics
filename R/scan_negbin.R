# Functions in this file:
#   scan_negbin
#   gen_negbin_counts
#   sim_negbin_statistic
#   negbin_mcsim
#   negbin_overdispersion
#   negbin_score_terms
#   poisson_score_terms
#   score_zone_sums
#   negbin_calculations
#   negbin_score
#   negbin_increasing_calculations
#   negbin_increasing_score
#   convolute_numerator
#   convolute_denominator


# Main function ----------------------------------------------------------------

#' Calculate the negative binomial scan statistic.
#' 
#' Calculate the expectation-based negative binomial scan statistic by supplying 
#' a \code{data.table} of observed counts and pre-computed distribution 
#' parameters for each location and time. A p-value for the observed scan
#' statistic can be obtained by Monte Carlo simulation.
#' 
#' @param table A \code{data.table} with columns \code{location, duration, mu,
#'    theta, count}. The \code{location} column should consist of integers that 
#'    are unique to each location. The \code{duration} column should also 
#'    consist of integers, starting at 1 for the most recent time period and 
#'    increasing in reverse chronological order. 
#'    
#'    A negative binomial distribution parametrized by \eqn{\mu} and 
#'    \eqn{\theta} (columns \code{mu} and \code{theta} respectively) has 
#'    expected value \eqn{\mu} and variance \eqn{\mu+\mu^2/\theta}. The 
#'    parameter \eqn{\theta} is referred to as the \code{size} in 
#'    \code{\link[stats]{NegBinomial}}, and \code{theta} in 
#'    \code{\link[MASS]{negative.binomial}}. 
#' @param zones A \code{set} of zones, each zone itself a 
#'    set containing one or more locations of those found in \code{table}.
#' @param n_mcsim A non-negative integer; the number of replicate scan 
#'    statistics to generate in order to calculate a p-value.
#' @param version Which version of the negative binomial score scan statistic to 
#'    calculate: either "ordinary" (default) or "increasing". See details.
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
#'     \item{distribution}{The assumed distribution of the data; "negative 
#'                         binomial" in this case.}
#'     \item{type}{The type of scan statistic; "Expectation-based" in this 
#'                 case.}
#'     \item{zones}{The set of zones that was passed to the function as input.}
#'     \item{n_locations}{The number of locations in the data.}
#'     \item{n_zones}{The number of zones.}
#'     \item{max_duration}{The maximum outbreak/event/anomaly duration 
#'                         considered.}
#'    }
#' @export
#' @concept negative binomial negbin nbinom scanstatistic
#' @details For the expectation-based negative binomial scan statistic (Tango
#'    et al., 2011), the null hypothesis of no anomaly holds that the count 
#'    observed at each location \eqn{i} and duration \eqn{t} (the number of time 
#'    periods before present) has a negative binomial distribution with expected 
#'    value \eqn{\mu_{it}} and dispersion parameter \eqn{\theta_{it}}:
#'    \deqn{
#'      H_0 : Y_{it} \sim \textrm{NegBin}(\mu_{it}, \theta_{it}).
#'    }
#'    This holds for all locations \eqn{i = 1, \ldots, m} and all durations 
#'    \eqn{t = 1, \ldots,T}, with \eqn{T} being the maximum duration considered.
#'    The alternative hypothesis depends on the version used: if \code{version
#'    == "ordinary"}, then the alternative hypothesis states that there is a 
#'    space-time window \eqn{W} consisting of a spatial zone \eqn{Z \subset \{1, 
#'    \ldots, m\}} and a time window \eqn{D \subseteq \{1, \ldots, T\}} such 
#'    that the counts in this window have their expected values inflated by a 
#'    factor \eqn{q_W > 1} compared to the null hypothesis:
#'    \deqn{
#'    H_1 : Y_{it} \sim \textrm{NegBin}(q_W \mu_{it}, \theta_{it}), 
#'          ~~(i,t) \in W.
#'    }
#'    If \code{version == "increasing"}, \eqn{q_W} is instead increasing over
#'    time (decreasing with \code{duration}).
#'    For locations and durations outside of this window, counts are assumed to
#'    be distributed as under the null hypothesis. The sets \eqn{Z} considered 
#'    are those specified in the argument \code{zones}, while the maximum 
#'    duration \eqn{T} is taken as the maximum value in the column 
#'    \code{duration} of the input \code{table}. For each space-time window
#'    \eqn{W} considered, a score statistic is computed using the score function
#'    and Fisher information under the null hypothesis of no anomaly.
#'    The scan statistic is calculated as the maximum of these quantities over 
#'    all space-time windows. Point estimates of the parameters \eqn{\mu_{it}} 
#'    and \eqn{\theta_{it}} must be specified in the column \code{mu} and 
#'    \code{theta} of the argument \code{table} before this function is called.
#' @references 
#'    Tango, T., Takahashi, K. & Kohriyama, K. (2011), \emph{A space-time scan 
#'    statistic for detecting emerging outbreaks}, Biometrics 67(1), 106–115.
#' @examples
#' # Simple example
#' set.seed(1)
#' table <- scanstatistics:::create_table(list(location = 1:4, duration = 1:4),
#'                                         keys = c("location", "duration"))
#' table[, mu := 3 * location]
#' table[, theta := 2]
#' table[, count := rnbinom(.N, mu = mu, size = theta)]
#' table[location %in% c(1, 4) & duration < 3, 
#'       count :=  rnbinom(.N, mu = 2 * mu, size = theta)]
#' zones <- scanstatistics:::powerset_zones(4)
#' result1 <- scan_negbin(table, zones, 100, "ordinary")
#' result2 <- scan_negbin(table, zones, 100, "increasing")
scan_negbin <- function(table, zones, n_mcsim = 0, version = "ordinary") {
  validate_scan(table, 
                zones, 
                c("count", "mu", "duration", "location", "theta"))
  details <- list(table = table,
                            zones = zones, 
                            distribution = "negative binomial",
                            type = "expectation-based",
                            version = version)
  if (version == "increasing") {
    scanstatistic_object(negbin_increasing_calculations(table, zones), 
                         negbin_mcsim(table, zones, n_mcsim, version),
                         details)
  } else {
    scanstatistic_object(negbin_calculations(table, zones), 
                         negbin_mcsim(table, zones, n_mcsim, version),
                         details)
  }
}

# Simulation and hypothesis testing functions ----------------------------------

#' Randomly generate and add negative binomial counts to a table.
#' 
#' This function randomly generates counts from a negative binomial distribution
#' according to the parameters on each row of the input \code{data.table}, and 
#' adds the counts to a new column \code{count}. 
#' @param table A \code{data.table} with at least the columns \code{mu} and
#'    \code{theta}. The parameter \eqn{\theta} (theta) is the same as 
#'    \code{size} in \code{\link[stats]{rnbinom}}.
#' @return The same table, with a new column \code{count}.
#' @importFrom stats rnbinom 
#' @importFrom stats rpois
#' @keywords internal
gen_negbin_counts <- function(table) {
  table[is.finite(theta), 
        count := as.integer(rnbinom(.N, mu = mu, size = theta))]
  table[is.infinite(theta), count := rpois(.N, mu)][]
}

#' Simulate a single negative binomial score scan statistic.
#' 
#' Simulate negative binomial-distributed data according to the supplied 
#' parameters and calculate the value of the score scan statistic, according to 
#' the specified model.
#' @inheritParams gen_negbin_counts
#' @inheritParams partition_zones
#' @param wstat_fun The function that calculates the statistic for each window.
#' @return A scalar; the scan statistic for the simulated data.
#' @importFrom magrittr %>%
#' @keywords internal
sim_negbin_statistic <- function(table, zones, wstat_fun) {
  table[, .(mu, theta), by = .(location, duration)] %>%
    gen_negbin_counts %>%
    negbin_overdispersion %>%
    wstat_fun(zones) %>%
    extract_scanstatistic
}

#' Monte Carlo simulation of negative binomial score scan statistics.
#' 
#' This function generates \code{n_mcsim} negative binomial-distributed 
#' data sets according to the parameters in the input table, and calculates the 
#' value of the score scan statistic for each generated data set using the 
#' supplied \code{zones}. The score can be calculated either according to the 
#' ordinary cluster model or the increasing outbreak/event/anomaly model.
#' @inheritParams gen_negbin_counts
#' @inheritParams partition_zones
#' @param n_mcsim A positive integer; the number of replicate scan 
#'    statistics to generate.
#' @param version Either "ordinary" (default) or "increasing".
#' @return A numeric vector of length \code{n_mcsim}.
#' @importFrom magrittr %>%
#' @keywords internal
negbin_mcsim <- function(table, zones, n_mcsim, version = "ordinary") {
  if (version == "increasing") {
    window_stats <- negbin_increasing_calculations
  } else {
    window_stats <- negbin_calculations
  }
  replicate(n_mcsim,
            table[, .(mu, theta), by = .(location, duration)] %>%
              gen_negbin_counts %>%
              negbin_overdispersion %>%
              window_stats(zones) %>%
              extract_scanstatistic)
}

#' Computes the overdispersion parameter for a fitted negative binomial model.
#' 
#' Computes the overdispersion parameter \eqn{w=1+\mu/\theta} for a negative
#' binomial distribution parametrized by its expected value \eqn{\mu} and with 
#' variance \eqn{\mu+\mu^2/\theta}. The overdispersion is added as a new column 
#' to the input \code{data.table}, meaning that this function \code{modifies} 
#' its input.
#' @param table A \code{data.table} with columns \code{mu, theta} and possibly
#'    others.
#' @return The same table, with a new column \code{overdispersion}.
#' @keywords internal
negbin_overdispersion <- function(table) {
  table[, overdispersion := 1 + mu / theta][]
}


#' Computes the numerator and denominator terms for the NegBin score.
#' 
#' This function calculates the terms found in the numerator and denominator 
#' sums for the ordinary version of the negative binomial scan statistic.
#' @param table A \code{data.table} with columns \code{location, duration, mu,
#'    overdispersion, count}.
#' @return A \code{data.table} with columns \code{location, duration, num, 
#'    denom}.
#' @keywords internal
negbin_score_terms <- function(table) {
  table[, 
        .(num = sum((count - mu) / overdispersion),
          denom = sum(mu / overdispersion)),
        by = .(location, duration)]
}

#' Computes the numerator and denominator terms for the Poisson score.
#' 
#' This function calculates the terms found in the numerator and denominator 
#' sums for the ordinary version of the Poisson scan statistic (Tango et al. 
#' version). 
#' @param table A \code{data.table} with columns \code{location, duration, mu,
#'    count}.
#' @return A \code{data.table} with columns \code{location, duration, num, 
#'    denom}.
#' @keywords internal
poisson_score_terms <- function(table) {
  table[,
        .(num = sum((count - mu)),
          denom = sum(mu)),
        by = .(location, duration)]
}

#' Sums the numerator and denominator terms over all locations in each zone.
#' 
#' Computes the sum of the numerator and denominator terms over all locations in
#' each zone, as part of the score calculation.
#' @param table A \code{data.table} with columns \code{location, duration, num,
#'    denom}; the output from \code{\link{negbin_score_terms}}.
#' @inheritParams partition_zones
#' @return A \code{data.table} with columns \code{zone, duration, num, denom}.
#' @importFrom magrittr %>%
#' @keywords internal
score_zone_sums <- function(table, zones) {
  table %>% 
    join_zones(zones = zones, keys = c("zone", "duration")) %>%
    zone_sum(sumcols = c("num", "denom"))
}


### Functions for ordinary model -----------------------------------------------

#' Calculate the ordinary NegBin score for each space-time window.
#' 
#' Calculate the ordinary negative binomial score for each space-time window, 
#' given the initial data of counts, expected values, and overdispersion 
#' parameters.
#' @inheritParams scan_negbin
#' @inheritParams partition_zones
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#' @importFrom magrittr %>%
#' @keywords internal
negbin_calculations <- function(table, zones) {
  table[, .(location, duration, count, mu, theta)] %>% 
    negbin_overdispersion %>%
    negbin_score_terms %>%
    score_zone_sums(zones) %>%
    negbin_score
}

#' Computes the ordinary score for each space-time window.
#' 
#' Computes the score statistic for each space-time window, assuming an 
#' ordinary outbreak/event/anomaly model and either a Poisson or a negative 
#' binomial distribution for the counts.
#' @param table A \code{data.table} with columns \code{zone, duration, num, 
#'    denom}.
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#' @keywords internal
negbin_score <- function(table) {
  table[,
        .(duration = duration,
          statistic = cumsum(num) / sqrt(cumsum(denom))),
        by = .(zone)]
}

### Functions for increasing model ---------------------------------------------

#' Calculate the increasing score for each space-time window.
#' 
#' Calculate the increasing score for each space-time window, given the initial 
#' data of counts, expected values, and overdispersion parameters.
#' @inheritParams scan_negbin
#' @inheritParams partition_zones
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#' @importFrom magrittr %>%
#' @keywords internal
negbin_increasing_calculations <- function(table, zones) {
  table[, .(location, duration, count, mu, theta)] %>% 
    negbin_overdispersion %>%
    negbin_score_terms %>%
    score_zone_sums(zones) %>%
    negbin_increasing_score
}

#' Calculate the increasing score for each space-time window.
#' 
#' Computes the score statistic for each space-time window, assuming the 
#' increasing type of outbreak/event/anomaly and either a Poisson or a negative 
#' binomial distribution for the counts.
#' @inheritParams negbin_score
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#' @keywords internal
negbin_increasing_score <- function(table) {
  table[,
    .(duration = duration,
      statistic = convolute_numerator(num, duration)
      / sqrt(convolute_denominator(denom, duration))), 
    by = .(zone)]
}

#' Computes the sum in the increasing score numerator.
#' 
#' @param x A vector of normalized counts summed over a single zone.
#' @param d A vector of outbreak durations considered.
#' @return A vector of length \code{length(d)}.
#' @keywords internal
convolute_numerator <- Vectorize(
  function(x, d) sum(d:1 * x[1:d]), vectorize.args = "d")

#' Computes the sum in the increasing score denominator.
#' @inheritParams convolute_numerator
#' @return A vector of length \code{length(d)}.
#' @keywords internal
convolute_denominator <- Vectorize(
  function(x, d) sum((d:1)^2 * x[1:d]), vectorize.args = "d")