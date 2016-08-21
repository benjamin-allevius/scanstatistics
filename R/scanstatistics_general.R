# Functions in this file:
#   extract_scanstatistic
#   extract_mlc
#   mc_pvalue
#   scanstatistic_object
#   print.scanstatistic
#   score_locations


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

#' Print a scanstatistic object.
#' 
#' Prints a scanstatistic object and returns it invisibly.
#' @param x A an object of class \code{scanstatistic}.
#' @param ... Further arguments passed to or from other methods.
#' @export
#' @keywords internal
print.scanstatistic <- function(x, ...) {
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
  invisible(x)
}

#' Score each location over zones and duration.
#' 
#' For each location, compute the average of the statistic calculated for each
#' space-time window that the location is included in, i.e. average the 
#' statistic over both zones and the maximum duration.
#' @param x An object of class \code{scanstatistic}.
#' @return A \code{data.table} with the following columns:
#'    \describe{
#'      \item{location}{The locations (as integers).}
#'      \item{total_score}{For each location, the sum of all window statistics 
#'                         that the location appears in.}
#'      \item{n_zones}{The number of spatial zones that the location appears 
#'                     in.}
#'      \item{score}{The total score divided by the number of zones and the 
#'                   maximum duration.}
#'      \item{relative_score}{The score divided by the maximum score.}
#' }
#' @export
#' @examples
#' # Simple example
#' set.seed(1)
#' table <- scanstatistics:::create_table(list(location = 1:4, duration = 1:4), 
#'                                         keys = c("location", "duration"))
#' table[, mean := 3 * location]
#' table[, count := rpois(.N, mean)]
#' table[location %in% c(1, 4) & duration < 3, count := rpois(.N, 2 * mean)]
#' zones <- scanstatistics:::powerset_zones(4)
#' result <- scan_poisson(table, zones, 100)
#' score_locations(result)
score_locations <- function(x) {
  tab <- data.table(location = seq_len(x$n_locations),
                    total_score = 0,
                    n_zones = 0)
  zone_scores <- x$observed[, .(score = sum(statistic)), by = zone]
  i <- 1
  for (z in x$zones) {
    locs <- unlist(z)
    tab[locs, total_score := total_score + zone_scores[locs, sum(score)]]
    tab[locs, n_zones := n_zones + 1]
    i <- i + 1
  }
  tab[, score := total_score / (n_zones * x$max_duration)]
  tab[, relative_score := score / max(score)]
  tab
}