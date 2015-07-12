

poisson_scanstatistic <- function(table, regions, n_replicates) {
  # input validation
  
  # Calculate statistics for observed data
  # extract max value
  observed_statistics <- poisson_calculations(table, regions)
  scan_obs <- extract_scanstatistic(observed_statistics)
  
  replicate_scanstats <- poisson_mcsim(table, regions, n_replicates)
  pval <- (1 + sum(replicate_scanstats > scan_obs)) / (1 + n_replicates)
  
  list(data = table,
       regions = regions,
       n_replicates = n_replicates,
       observed = observed_statistics,
       mlc = extract_mlc(observed_statistics),
       pvalue = pval)
}

#' Monte Carlo simulation of expectation-based Poisson scan statistics.
#' 
#' This function generates \code{n_replicates} Poisson-distributed data sets 
#' according to the parameters in the input table, and calculates the value of
#' the scan statistic for each generated data set using the supplied 
#' \code{regions}.
#' @param table A \code{data.table} with columns \code{location, duration, 
#'    mean}.
#' @inheritParams partition_regions
#' @param n_replicates A positive integer; the number of replicate scan 
#'    statistics to generate.
#' @return A numeric vector of length \code{n_replicates}.
#' @importFrom magrittr %>%
poisson_mcsim <- function(table, regions, n_replicates = 999L) {
  foreach::foreach(i = seq(n_replicates), 
                          .combine = c, 
                          .inorder = FALSE) %do% {
    table[, .(mean), by = .(location, duration)] %>%
      generate_poisson_counts %>% 
      poisson_calculations(regions = regions) %>%
      extract_scanstatistic
  }
}

#' Randomly generate and add Poisson counts to a table.
#' 
#' This function randomly generates counts from a Poisson distribution according 
#' to the parameters on each row of the input \code{data.table}, and adds the 
#' counts to a new column \code{count}. 
#' @param table A \code{data.table} with at least the column \code{mean}, 
#'    corresponding to the parameter \code{lambda} in 
#'    \code{\link[stats]{rpois}}.
generate_poisson_counts <- function(table) {
  table[, count := rpois(.N, lambda = mean)][]
}

#' Calculate the expectation-based Poisson statistic for each space-time window.
#' 
#' Calculate the expectation-based Poisson statistic for each space-time window,
#' given the initial data of counts and means.
#' @param table A \code{data.table} with columns \code{location, duration, 
#'    count, mean}.
#' @param regions A \code{list} or \code{set} of regions, each region itself a 
#'    set containing one or more locations of those found in \code{table}.
#' @return A \code{data.table} with columns \code{region, duration, statistic}.
#'    The column \code{statistic} contains the logarithm of the statistic for 
#'    each region-duration combination.
#' @importFrom magrittr %>%
poisson_calculations <- function(table, regions) {
  table %>% 
    region_joiner(regions = regions, keys = c("region", "duration")) %>%
    region_sum(sumcols = c("count", "mean")) %>%
    cumsum_duration(sumcols = c("count", "mean"), bycols = c("region")) %>%
    poisson_relrisk %>%
    poisson_statistic
}

#' Calculate the expectation-based Poisson statistic for each space-time window.
#' 
#' This function calculates the logarithm of the expectation-based Poisson 
#' statistic for each space-time window, given already calculated relative risks
#' and aggregate counts and means.
#' @param table A \code{data.table} with columns \code{region, duration, count,
#'    mean, relrisk}. The columns \code{count} and \code{mean} contain the sums
#'    of the counts and means for the locations inside each region and up to the
#'    given duration. The column \code{relrisk} contains the maximum likelihood
#'    estimate for the relative risk for each region-duration combination.
#' @return A \code{data.table} with columns \code{region, duration, statistic}.
#'    The column \code{statistic} contains the logarithm of the statistic for 
#'    each region-duration combination.
poisson_statistic <- function(table) {
  table[, .(statistic = log(relrisk) * count - (relrisk - 1) * mean),
        by = .(region, duration)]
}

#' Add column for relative risk MLE to table with aggregate counts and means.
#' 
#' This function adds a column for the maximum likelihood estimate for the
#' relative risk (assumed to be 1 or greater). The table should contain the
#' aggregates (sums) of the counts and means for each region and duration.
#' @param table A \code{data.table} with columns \code{region, duration, count,
#'    mean}. The latter two are the sums of counts and means for the locations
#'    comprising the region and up to the given duration.
#' @return The same table, with an extra column \code{relrisk}.
poisson_relrisk <- function(table) {
  table[, relrisk := max(1, count / mean), by = .(region, duration)]
}