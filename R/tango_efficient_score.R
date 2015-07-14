### General functions ----------------------------------------------------------

nbinom_score_scanstatistic <- function(
  table, regions, n_replicates, type = "hotspot") {
  # input validation
  
  if (type == "outbreak") {
    window_stats <- outbreak_calculations
  } else {
    window_stats <- hotspot_calculations
  }
  
  # Calculate statistics for observed data
  # extract max value
  observed_statistics <- window_stats(table, regions)
  scan_obs <- extract_scanstatistic(observed_statistics)
  
  replicate_scanstats <- nbinom_mcsim(table, regions, n_replicates, type)
  pval <- (1 + sum(replicate_scanstats > scan_obs)) / (1 + n_replicates)
  
  list(data = table,
       regions = regions,
       n_replicates = n_replicates,
       replicates = replicate_scanstats,
       observed = observed_statistics,
       mlc = extract_mlc(observed_statistics),
       pvalue = pval,
       replicates = replicate_scanstats)
}

#' Monte Carlo simulation of negative binomial efficient score scan statistics.
#' 
#' This function generates \code{n_replicates} negative binomial-distributed 
#' data sets according to the parameters in the input table, and calculates the 
#' value of the efficient score scan statistic for each generated data set using 
#' the supplied \code{regions}. The score can be calculated either according to
#' the hotspot cluster model or the emerging outbreak model.
#' @inheritParams generate_nbinom_counts
#' @inheritParams partition_regions
#' @param n_replicates A positive integer; the number of replicate scan 
#'    statistics to generate.
#' @return A numeric vector of length \code{n_replicates}.
#' @importFrom magrittr %>%
nbinom_mcsim <- function(table, regions, n_replicates, type = "hotspot") {
  if (type == "outbreak") {
    window_stats <- outbreak_calculations
  } else {
    window_stats <- hotspot_calculations
  }
  foreach::foreach(i = seq(n_replicates), 
                   .combine = c, 
                   .inorder = FALSE) %do% {
    table[, .(mean, phi), by = .(location, duration)] %>%
      generate_nbinom_counts %>%
      compute_nbinom_overdispersion %>%
      window_stats(regions) %>%
      extract_scanstatistic
  }
}

#' Randomly generate and add negative binomial counts to a table.
#' 
#' This function randomly generates counts from a negative binomial distribution
#' according to the parameters on each row of the input \code{data.table},
#' and adds the counts to a new column \code{count}. 
#' @param table A \code{data.table} with at least the columns \code{mean} and
#'    \code{phi}. The parameter \eqn{\phi} (phi) is the same as \code{size} in
#'    \code{\link[stats]{rnbinom}}.
#' @return The same table, with a new column \code{count}.
generate_nbinom_counts <- function(table) {
  table[, count := rnbinom(.N, mu = mean, size = phi)][]
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
compute_nbinom_overdispersion <- function(table) {
  table[, overdispersion := 1 + mean / phi][]
}


#' Computes the numerator and denominator terms for the hotspot efficient score.
#' 
#' This function calculates the terms found in the numerator and denominator 
#' sums for the hotspot version of the negative binomial efficient score. 
#' @param table A \code{data.table} with columns \code{location, duration, mean,
#'    overdispersion, count}. If \eqn{\mu} is the mean of the negative binomial 
#'    distribution and \eqn{\phi} is the parameter such that the variance of the 
#'    distribution is \eqn{\mu+\mu^2/\phi}, the overdispersion is given by 
#'    \eqn{1+\mu/\phi}. The parameter \eqn{\phi} is referred to as the 
#'    \code{size} in \code{\link[stats]{NegBinomial}}, and \code{theta} 
#'    in \code{\link[MASS]{negative.binomial}}.
#' @return A \code{data.table} with columns \code{location, duration, num, 
#'    denom}.
efficient_score_terms_nbinom <- function(table) {
  table[, 
        .(num = sum((count - mean) / overdispersion),
          denom = sum(mean / overdispersion)),
        by = .(location, duration)]
}

#' Computes the numerator and denominator terms for the hotspot efficient score.
#' 
#' This function calculates the terms found in the numerator and denominator 
#' sums for the hotspot version of the Poisson efficient score. 
#' @param table A \code{data.table} with columns \code{location, duration, mean,
#'    count}.
#' @return A \code{data.table} with columns \code{location, duration, num, 
#'    denom}.
efficient_score_terms_poisson <- function(table) {
  table[,
        .(num = sum((count - mean)),
          denom = sum(mean)),
        by = .(location, duration)]
}

#' Sums the numerator and denominator terms over all locations in each region.
#' 
#' Computes the sum of the numerator and denominator terms over all locations in
#' each region, as part of the efficient score calculation.
#' @param table A \code{data.table} with columns \code{location, duration, num,
#'    denom}; the output from \code{\link{efficient_score_terms_nbinom}}.
#' @inheritParams partition_regions
#' @return A \code{data.table} with columns \code{region, duration, num, denom}.
#' @importFrom magrittr %>%
efficient_score_region_sums <- function(table, regions) {
  table %>% 
    region_joiner(regions = regions, keys = c("region", "duration")) %>%
    region_sum(sumcols = c("num", "denom"))
}


### Functions for hotspot model ------------------------------------------------

#' Calculate the hotspot efficient score for each space-time window.
#' 
#' Calculate the hotspot efficient score for each space-time window, given the 
#' initial data of counts, means, and overdispersion parameters.
#' @inheritParams efficient_score_terms_nbinom
#' @inheritParams partition_regions
#' @return A \code{data.table} with columns \code{region, duration, statistic}.
#' @importFrom magrittr %>%
hotspot_calculations <- function(table, regions) {
  table %>% 
    efficient_score_terms_nbinom %>%
    efficient_score_region_sums(regions) %>%
    hotspot_efficient_score
}

#' Computes the hotspot efficient score for each space-time window.
#' 
#' Computes the efficient score statistic for each space-time window, assuming a 
#' hotspot outbreak model and either a Poisson or a negative binomial 
#' distribution for the counts.
#' @param table A \code{data.table} with columns \code{region, duration, num, 
#'    denom}.
#' @return A \code{data.table} with columns \code{region, duration, statistic}.
hotspot_efficient_score <- function(table) {
  table[,
        .(duration = duration,
          statistic = cumsum(num) / sqrt(cumsum(denom))),
        by = .(region)]
}

### Functions for outbreak model -----------------------------------------------

#' Calculate the outbreak efficient score for each space-time window.
#' 
#' Calculate the outbreak efficient score for each space-time window, given the 
#' initial data of counts, means, and overdispersion parameters.
#' @inheritParams efficient_score_terms_nbinom
#' @inheritParams partition_regions
#' @return A \code{data.table} with columns \code{region, duration, statistic}.
#' @importFrom magrittr %>%
outbreak_calculations <- function(table, regions) {
  table %>% 
    efficient_score_terms_nbinom %>%
    efficient_score_region_sums(regions) %>%
    outbreak_efficient_score
}

#' Calculate the outbreak efficient score for each space-time window.
#' 
#' Computes the efficient score statistic for each space-time window, assuming 
#' an (emergent) outbreak model and either a Poisson or a negative binomial 
#' distribution for the counts.
#' @inheritParams hotspot_efficient_score
#' @return A \code{data.table} with columns \code{region, duration, statistic}.
outbreak_efficient_score <- function(table) {
  table[,
    .(duration = duration,
      statistic = convolute_numerator(num, duration)
      / sqrt(convolute_denominator(denom, duration))), 
    by = .(region)]
}

#' Computes the sum in the outbreak efficient score numerator.
#' 
#' @param x A vector of normalized counts summed over a single region.
#' @param d A vector of outbreak durations considered.
#' @return A vector of length \code{length(d)}.
convolute_numerator <- Vectorize(
  function(x, d) sum(d:1 * x[1:d]), vectorize.args = "d")

#' Computes the sum in the outbreak efficient score denominator.
#' @inheritParams convolute_numerator
#' @return A vector of length \code{length(d)}.
convolute_denominator <- Vectorize(
  function(x, d) sum((d:1)^2 * x[1:d]), vectorize.args = "d")