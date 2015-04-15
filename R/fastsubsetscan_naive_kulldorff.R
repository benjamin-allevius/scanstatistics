
#' Score and minimal stream subset for the Naive Kulldorff method.
#' 
#' Calculates the score for each combination of region and duration, and the 
#' minimal subset of locations that contribute to the score, according to the
#' Naive Kulldorff method. See Neill et. al. (2013) section 3.2.
#' @param scores A \code{data.table} with columns \code{region, duration, stream, 
#'    score}.
#' @return A \code{data.table} with columns \code{region, duration, stream, 
#'    score, included_streams}. The colum \code{score} contain the sum of the 
#'    input score over all data streams, for each region and duration. The 
#'    column \code{included_streams} contain those data streams that made a 
#'    positive contribution to this sum.
score_minimal_stream_subset <- function(scores) {
  scores[score > 0, 
         .(score = sum(score), included_streams = list(stream)), 
         by = .(region, duration)]
}

#' Calculates the multivariate scan statistic by the naive Kulldorff method.
#' 
#' Calculates the multivariate scan statistic by the naive Kulldorff method, for 
#' all stream-region-duration combinations, according to the distribution used 
#' in the supplied functions.
#' @param counts A \code{data.table} with columns \code{region, duration, 
#'    stream, count, baseline}, and possibly others.
#' @param regions A \code{list} or \code{set} of regions, 
#'    each region itself a set containing one or more locations of those found 
#'    in \code{counts}.
#' @inheritParams region_apply
#' @param aggregate_CB A function for calculating the aggregate counts and 
#'    baselines over all event durations for a given expectation-based scan 
#'    statistic.
#' @inheritParams score_EB
#' @return A \code{data.table} with columns \code{region, duration, score,
#'    included_streams}.
#' @importFrom magrittr %>%
naive_kulldorff_general <- function(counts, 
                                    regions, 
                                    region_partition,
                                    aggregate_CB, 
                                    score_function) {
  counts %>%
    aggregate_CB %>%
    region_apply(region_partition = region_partition,
                 f = aggregate_per_stream,
                 keys = c("region", "duration", "stream")) %>%
    score_EB(score_function = score_function) %>%
    score_minimal_stream_subset
}

# Expectation-based Poisson ----------------------------------------------------

#' Calculates the EBP multivariate scan statistic by the naive Kulldorff method.
#' 
#' This function calculates the multivariate scan statistic by the naive 
#' Kulldorff method, using the expectation-based Poisson score function.
#' @inheritParams naive_kulldorff_general
naive_kulldorff_poisson <- function(counts, regions, region_partition) {
  naive_kulldorff_general(
    counts, regions, region_partition, aggregate_CB_poisson, score_fun_EBP)
}


# Expectation-based Gaussian ---------------------------------------------------

#' Calculates the EBG multivariate scan statistic by the naive Kulldorff method.
#' 
#' This function calculates the multivariate scan statistic by the naive 
#' Kulldorff method, using the expectation-based Gaussian score function.
#' @inheritParams naive_kulldorff_general
naive_kulldorff_gaussian <- function(counts, regions, region_partition) {
  naive_kulldorff_general(
    counts, regions, region_partition, aggregate_CB_gaussian, score_fun_EBG)
}


# Expectation-based Exponential ------------------------------------------------

#' Calculates the EBE multivariate scan statistic by the naive Kulldorff method.
#' 
#' This function calculates the multivariate scan statistic by the naive 
#' Kulldorff method, using the expectation-based exponential score function.
#' @inheritParams naive_kulldorff_general
naive_kulldorff_exponential <- function(counts, regions, region_partition) {
  naive_kulldorff_general(
    counts, regions, region_partition, aggregate_CB_exponential, score_fun_EBE)
}