
#' Calculates the multivariate scan statistic by the naive Kulldorff method.
#' 
#' Calculates the multivariate scan statistic by the naive Kulldorff method, for 
#' all stream-region-duration combinations, according to the given distribution.
#' @param counts A \code{data.table} with columns \code{region, duration, 
#'    stream, count, baseline}, and possibly others.
#' @param distribution The distribution to be used; one of "poisson", "gaussian"
#'    (or "normal"), or "exponential".
#' @param regions A \code{list} or \code{set} of regions, 
#'    each region itself a set containing one or more locations of those found 
#'    in \code{counts}.
#' @inheritParams MBSS
#' @return A \code{data.table} with columns \code{region, duration, score,
#'    included_streams}.
#' @importFrom magrittr %>%
naive_kulldorff <- function(counts, 
                            distribution,
                            regions, 
                            n_partitions = "auto") {
  if (n_partitions == "auto") {
    n_partitions <- auto_region_partition_size(regions)
  }
  region_partition <- partition_regions(regions, n_parts = n_partitions)
  
  # Define aggregation and score functions based on distribution
  initial_aggregation_fun <- dispatch_initial_aggregation(distribution)
  score_fun <- dispatch_score_function(distribution)
  
  counts %>%
    initial_aggregation_fun %>%
    region_apply(region_partition = region_partition,
                 f = aggregate_per_stream,
                 keys = c("region", "duration", "stream")) %>%
    expectation_based_score(score_function = score_fun,
                            region_as_list = FALSE) %>%
    score_minimal_stream_subset(region_as_list = FALSE)
}