#' Calculate a scan statistic and its corresponding p-value, MLC, etc
#' @param scanstatistic The type of scan statistic. Can be "poisson", 
#'    "nb-hotspot", "nb-emerging", "zip". See details.
#' @param table A \code{data.table} with columns \code{duration, location, 
#'    count} and columns for distribution parameters that match the chosen type
#'    of scan statistic. See details.
#' @param zones A \code{set} of \code{set}s, the inner sets containing integer
#'    values that match those in the column \code{location} of the \code{table}
#'    argument of this function.
#' @param n_replicates The number of Monte Carlo replications to perform in 
#'    order to calculate a p-value.
#' @param ... Arguments passed to 
#' @export
calculate_scanstatistic <- function(scanstatistic = "poisson",
                                    table, 
                                    zones, 
                                    n_replicates, ...) {
  # Input validation
  allowed_types <- c("poisson", "nb-hotspot", "nb-emerging", "zip")
  if (scanstatistic %notin% allowed_types) {
    stop(paste0("Argument scanstatistic must be one of ",
                toString(allowed_types), "."))
  }
  
  # Calculate the chosen type of scanstatistic
  if (scanstatistic == "poisson") {
    observed_statistics <- poisson_calculations(table, zones)
    replicate_scanstats <- poisson_mcsim(table, zones, n_replicates)
  } else if (scanstatistic == "nb-hotspot") {
      observed_statistics <- nb_hotspot_calculations(table, zones)
      replicate_scanstats <- nb_mcsim(table, zones, n_replicates, "hotspot")
  } else if (scanstatistic == "nb-emerging") {
      observed_statistics <- nb_emerging_calculations(table, zones)
      replicate_scanstats <- nb_mcsim(table, zones, n_replicates, "emerging")
  } else if (scanstatistic == "zip") {
    maxdur <- table[, max(duration)]
    observed_statistics <- zip_calculations(table, zones, maxdur = maxdur, ...)
    replicate_scanstats <- zip_mcsim(
      table, zones, n_replicates, maxdur = maxdur, ...)
  }
  
  scan_obs <- extract_scanstatistic(observed_statistics)
  pval <- mc_pvalue(scan_obs, replicate_scanstats)
  mlc <- extract_mlc(observed_statistics)
  
  list(observed_scanstatistic = mlc[, statistic],
       pvalue = pval,
       outbreak_zone = mlc[, zone],
       outbreak_locations = get_zone(mlc[, zone], zones),
       outbreak_duration = mlc[, duration],
       observed_statistics = observed_statistics,
       replicate_scanstatistics = replicate_scanstats)
}


#' Extract value of scan statistic from per-window statistics.
#' 
#' This function extracts the value of the scan statistic, which is the maximum
#' of the statistics calculated for each spatial or space-time window.
#' @param table A \code{data.table} with a column \code{statistic}, which should
#'    correspond to the statistic calculated for each spatial or space-time 
#'    window (given as other columns).
#' @return The maximum value of the column \code{statistic}.
#' @keywords internal
extract_scanstatistic <- function(table) {
  table[, max(statistic)]
}

#' Extract the most likely cluster (MLC) and the value of the scans statistic.
#' 
#' This function extracts the most likely cluster, which is the spatial or 
#' spatiotemporal window that corresponds to the scan statistic. It also returns
#' the value of the scan statistc.
#' @inheritParams extract_scanstatistic
#' @return The row of the input table with the highest value of the column 
#'    \code{statistic}.
#' @keywords internal
extract_mlc <- function(table) {
  table[which.max(statistic), ]
}

#' Calculate the Monte Carlo p-value for a scan statistic.
#' @param observed A scalar; the observed value of the scan statistic.
#' @param replicates A vector of Monte Carlo replicates of the scan statistic.
#' @return A scalar; the p-value corresponding to the observed scan statistic.
#' @keywords internal
mc_pvalue <- function(observed, replicates) {
  (1 + sum(replicates > observed)) / (1 + length(replicates))
}

#' Creates an S3 object of class scanstatistic.
#' @keywords internal
scanstatistic_object <- function(observed, simulated, details) {
  statistic <- extract_scanstatistic(observed)
  pval <- mc_pvalue(statistic, replicated)
  mlc <- extract_mlc(observed)
  
  structure(list(observed = observed,
                 replicated = unlist(replicated),
                 mlc = mlc,
                 pvalue = pval,
                 distribution = details$distribution,
                 type = details$type,
                 n_locations = length(details$table[, unique(location)]),
                 n_zones = length(details$zones),
                 n_maxduration = details$table[, max(duration)],
                 class = "scanstatistic"))
}


print.scanstatistic <- function(x) {
  cat(paste0(
    "A scan statistic assuming a ", 
    x$distribution,
    " distribution for the data was run on a dataset consisting of ", 
    x$n_locations, 
    " locations, making up ",
    x$n_zones, 
    " zones. The maximum outbreak/event/anomaly duration considered was ",
    x$n_maxduration, "."))
}