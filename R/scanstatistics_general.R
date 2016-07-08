# Functions in this file:
#   extract_scanstatistic
#   extract_mlc
#   mc_pvalue
#   scanstatistic_object
#   print.scanstatistic


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
  if (length(replicates) == 0) {
    return(NULL)
  } else {
    return((1 + sum(replicates > observed)) / (1 + length(replicates)))
  }
}

#' Creates an S3 object of class scanstatistic.
#' @keywords internal
scanstatistic_object <- function(observed, simulated, details) {
  statistic <- extract_scanstatistic(observed)
  pval <- mc_pvalue(statistic, simulated)
  mlc <- extract_mlc(observed)
  
  structure(list(observed = observed,
                 replicated = unlist(simulated),
                 mlc = mlc,
                 pvalue = pval,
                 distribution = details$distribution,
                 type = details$type,
                 zones = details$zones,
                 n_locations = length(details$table[, unique(location)]),
                 n_zones = length(details$zones),
                 max_duration = details$table[, max(duration)]),
            class = "scanstatistic")
}


print.scanstatistic <- function(x) {
  cat(paste0(
    "Data distribution:                ", x$distribution, "\n",
    "Type of scan statistic:           ", x$type, "\n",
    "Number of locations considered:   ", x$n_locations, "\n",
    "Maximum duration considered:      ", x$max_duration, "\n",
    "Number of spatial zones:          ", x$n_zones, "\n",
    "p-value of observed statistic:    ", round(x$pvalue, 3), "\n",
    "Number of Monte Carlo replicates: ", length(x$replicated), "\n",
    "Most likely event duration:       ", x$mlc$duration, "\n",
    "ID of locations in most likely cluster: ", 
    toString(get_zone(x$mlc$zone, x$zones)))
    )
}