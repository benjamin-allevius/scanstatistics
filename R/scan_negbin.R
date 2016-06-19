# Main function ----------------------------------------------------------------

#' @param table A \code{data.table} with columns \code{location, duration, mean,
#'    overdispersion, count}. If \eqn{\mu} is the mean of the negative binomial 
#'    distribution and \eqn{\phi} is the parameter such that the variance of the 
#'    distribution is \eqn{\mu+\mu^2/\phi}, the overdispersion is given by 
#'    \eqn{1+\mu/\phi}. The parameter \eqn{\phi} is referred to as the 
#'    \code{size} in \code{\link[stats]{NegBinomial}}, and \code{theta} 
#'    in \code{\link[MASS]{negative.binomial}}.
#' @param n_replicates A positive integer; the number of replicate scan 
#'    statistics to generate. 
#' @param version Which version of the negative binomial score scan statistic to 
#'    calculate: either "ordinary" (default) or "increasing". See details.
#' @inheritParams partition_zones
#' @return An object of class \code{scanstatistics}.
scan_negbin <- function(table, zones, n_replicates = 0, version = "ordinary") {
  details <- list(table = table,
                            zones = zones, 
                            distribution = "negative binomial",
                            type = "expectation-based",
                            version = version)
  if (version == "increasing") {
    scanstatistic_object(negbin_increasing_calculations(table, zones), 
                         negbin_mcsim(table, zones, n_replicates, version),
                         details)
  } else {
    scanstatistic_object(negbin_calculations(table, zones), 
                         negbin_mcsim(table, zones, n_replicates, version),
                         details)
  }
}

# Simulation and hypothesis testing functions ----------------------------------

#' Randomly generate and add negative binomial counts to a table.
#' 
#' This function randomly generates counts from a negative binomial distribution
#' according to the parameters on each row of the input \code{data.table}, and 
#' adds the counts to a new column \code{count}. 
#' @param table A \code{data.table} with at least the columns \code{mean} and
#'    \code{phi}. The parameter \eqn{\phi} (phi) is the same as \code{size} in
#'    \code{\link[stats]{rnbinom}}.
#' @return The same table, with a new column \code{count}.
#' @importFrom stats rnbinom 
#' @importFrom stats rpois
#' @keywords internal
gen_negbin_counts <- function(table) {
  table[is.finite(phi), count := as.integer(rnbinom(.N, mu = mean, size = phi))]
  table[is.infinite(phi), count := rpois(.N, mean)][]
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
  table[, .(mean, phi), by = .(location, duration)] %>%
    gen_negbin_counts %>%
    negbin_overdispersion %>%
    wstat_fun(zones) %>%
    extract_scanstatistic
}

#' Monte Carlo simulation of negative binomial score scan statistics.
#' 
#' This function generates \code{n_replicates} negative binomial-distributed 
#' data sets according to the parameters in the input table, and calculates the 
#' value of the score scan statistic for each generated data set using the 
#' supplied \code{zones}. The score can be calculated either according to the 
#' ordinary cluster model or the increasing outbreak/event/anomaly model.
#' @inheritParams gen_negbin_counts
#' @inheritParams partition_zones
#' @param n_replicates A positive integer; the number of replicate scan 
#'    statistics to generate.
#' @param version Either "ordinary" (default) or "increasing".
#' @return A numeric vector of length \code{n_replicates}.
#' @importFrom magrittr %>%
#' @importFrom foreach foreach
#' @keywords internal
negbin_mcsim <- function(table, zones, n_replicates, version = "ordinary") {
  if (version == "increasing") {
    window_stats <- negbin_increasing_calculations
  } else {
    window_stats <- negbin_calculations
  }
  replicate(n_replicates,
            table[, .(mean, phi), by = .(location, duration)] %>%
              gen_negbin_counts %>%
              negbin_overdispersion %>%
              window_stats(zones) %>%
              extract_scanstatistic)
}

#' Computes the overdispersion parameter for a fitted negative binomial model.
#' 
#' Computes the overdispersion parameter \eqn{w=1+\mu/\phi} for a negative
#' binomial distribution parametrized by its mean \eqn{\mu} and with variance
#' \eqn{\mu+\mu^2/\phi}. The overdispersion is added as a new column to the 
#' input \code{data.table}, meaning that this function \code{modifies} its 
#' input.
#' @param table A \code{data.table} with columns \code{mean, phi} and possibly
#'    others.
#' @return The same table, with a new column \code{overdispersion}.
#' @keywords internal
negbin_overdispersion <- function(table) {
  table[, overdispersion := 1 + mean / phi][]
}


#' Computes the numerator and denominator terms for the NegBin score
#' 
#' This function calculates the terms found in the numerator and denominator 
#' sums for the ordinary version of the negative binomial scan statistic.
#' @inheritParams scan_negbin
#' @return A \code{data.table} with columns \code{location, duration, num, 
#'    denom}.
#' @keywords internal
negbin_score_terms <- function(table) {
  table[, 
        .(num = sum((count - mean) / overdispersion),
          denom = sum(mean / overdispersion)),
        by = .(location, duration)]
}

#' Computes the numerator and denominator terms for the Poisson score.
#' 
#' This function calculates the terms found in the numerator and denominator 
#' sums for the ordinary version of the Poisson scan statistic (Tango et al. 
#' version). 
#' @param table A \code{data.table} with columns \code{location, duration, mean,
#'    count}.
#' @return A \code{data.table} with columns \code{location, duration, num, 
#'    denom}.
#' @keywords internal
poisson_score_terms <- function(table) {
  table[,
        .(num = sum((count - mean)),
          denom = sum(mean)),
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
    zone_joiner(zones = zones, keys = c("zone", "duration")) %>%
    zone_sum(sumcols = c("num", "denom"))
}


### Functions for ordinary model ------------------------------------------------

#' Calculate the ordinary NegBin score for each space-time window.
#' 
#' Calculate the ordinary negative binomial score for each space-time window, 
#' given the initial data of counts, means, and overdispersion parameters.
#' @inheritParams scan_negbin
#' @inheritParams partition_zones
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#' @importFrom magrittr %>%
#' @keywords internal
negbin_calculations <- function(table, zones) {
  table %>% 
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

### Functions for increasing model -----------------------------------------------

#' Calculate the increasing score for each space-time window.
#' 
#' Calculate the increasing score for each space-time window, given the initial 
#' data of counts, means, and overdispersion parameters.
#' @inheritParams scan_negbin
#' @inheritParams partition_zones
#' @return A \code{data.table} with columns \code{zone, duration, statistic}.
#' @importFrom magrittr %>%
#' @keywords internal
negbin_increasing_calculations <- function(table, zones) {
  table %>% 
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