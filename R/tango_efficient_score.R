### General functions ----------------------------------------------------------

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
compute_nb_overdispersion <- function(table) {
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
efficient_score_terms_nbin <- function(table) {
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
#'    denom}; the output from \code{\link{efficient_score_terms_nbin}}.
#' @inheritParams partition_regions
#' @return A \code{data.table} with columns \code{region, duration, num, denom}.
#' @importFrom magrittr %>%
efficient_score_region_sums <- function(table, regions) {
  table %>% 
    region_joiner(regions = regions, keys = c("region", "duration")) %>%
    region_sum(sumcols = c("num", "denom"))
}


### Functions for hotspot model ------------------------------------------------

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