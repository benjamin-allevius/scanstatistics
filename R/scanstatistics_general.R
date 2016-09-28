# Functions in this file:
#   extract_scanstatistic
#   extract_mlc
#   mc_pvalue
#   scanstatistic_object
#   print.scanstatistic
#   score_locations
#   top_clusters
#   validate_colnames
#   validate_values
#   validate_zones
#   validate_scan


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

#' Calculate the Monte Carlo \eqn{p}-value for a scan statistic.
#' 
#' Given an observed scan statistic \eqn{\lambda^*} and a vector of replicate scan 
#' statistics \eqn{\lambda_i}, \eqn{i=1,\ldots,R}, calculate the Monte carlo 
#' \eqn{p}-value as
#' \deqn{
#'  \frac{1 + \sum_{i=1}^R \mathrm{I}(\lambda_i > \lambda^*)}{1 + R}
#' }
#' The function is vectorized, so multiple \eqn{p}-values can be calculated if
#' several scan statistics (e.g. statistics from secondary clusters) are 
#' supplied.
#' @param observed A scalar containing the observed value of the scan statistic,
#'    or a vector of observed values from secondary clusters.
#' @param replicates A vector of Monte Carlo replicates of the scan statistic.
#' @return The \eqn{p}-value or \eqn{p}-values corresponding to the observed 
#'    scan statistic(s).
#' @export
mc_pvalue <- function(observed, replicates) {
  if (length(replicates) == 0) {
    return(NULL)
  } else {
    f <- Vectorize(
      function(y) {
        (1 + sum(replicates > y)) / (1 + length(replicates))
        }
    )
    
    return(f(observed))
  }
}

#' Creates an S3 object of class scanstatistic.
#' @param observed A data table containing columns location, duration, 
#'    statistic, and possibly others.
#' @param simulated A numeric vector of replicated scan statistics.
#' @param details A list containing details about the data and scan statistic.
#' @keywords internal
scanstatistic_object <- function(observed, simulated, details) {
  statistic <- extract_scanstatistic(observed)
  pval <- mc_pvalue(statistic, simulated)
  mlc <- extract_mlc(observed)
  
  structure(list(observed = observed[order(-statistic), ],
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
    "Number of Monte Carlo replicates: ", length(x$replicated), "\n",
    "p-value of observed statistic:    ", ifelse(is.null(x$pvalue), 
                                                 "NULL",
                                                 round(x$pvalue, 3)), "\n",
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
#' table[, mu := 3 * location]
#' table[, count := rpois(.N, mu)]
#' table[location %in% c(1, 4) & duration < 3, count := rpois(.N, 2 * mu)]
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
    tab[z, total_score := total_score + zone_scores[z, sum(score)]]
    tab[z, n_zones := n_zones + 1]
    i <- i + 1
  }
  tab[, score := total_score / (n_zones * x$max_duration)]
  tab[, relative_score := score / max(score)]
  tab
}

#' Get the top (non-overlappig) clusters.
#' 
#' Get the top \eqn{k} space-time clusters according to the statistic calculated
#' for each cluster (the maximum being the scan statistic). The default is to 
#' return the spatially non-overlapping clusters, i.e. those that do not have 
#' any locations in common.
#' @param x An object of class scanstatistics.
#' @param k An integer, the number of clusters to return
#' @param overlapping Logical; should the top clusters be allowed to overlap in
#'    the spatial dimension? The default is \code{FALSE}.
#' @return A \code{data.table} with at most \eqn{k} rows, with columns 
#'    \code{zone, duration, statistic}. 
#' @export
#' @examples 
#' set.seed(1)
#' table <- scanstatistics:::create_table(list(location = 1:4, duration = 1:4), 
#'                                         keys = c("location", "duration"))
#' table[, mu := 3 * location]
#' table[, count := rpois(.N, mu)]
#' table[location %in% c(1, 4) & duration < 3, count := rpois(.N, 2 * mu)]
#' zones <- scanstatistics:::powerset_zones(4)
#' result <- scan_poisson(table, zones, 0)
#' top_clusters(result, k = 4, overlapping = FALSE)
top_clusters <- function(x, k = 5, overlapping = FALSE) {
  if (overlapping) {
    return(x$observed[seq_len(k), ])
  } else {
    row_idx <- integer(k)
    seen_locations <- integer(0)
    n_added <- 0L
    i <- 1L
    while (n_added < k && i < nrow(x$observed)) {
      zone <- x$observed[i, zone]
      if (length(intersect(seen_locations, x$zones[[zone]])) == 0) {
        n_added <- n_added + 1L
        seen_locations <- c(seen_locations, x$zones[[zone]])
        row_idx[n_added] <- i
      }
      i <- i + 1L
    }
    return(x$observed[row_idx[row_idx > 0], ])
  }
}


#' Check that the input table has the right columns. Raises error if not.
#' @param table A data.table passed to a scan statistic function.
#' @param col_names A character vector of column names required.
#' @keywords internal
validate_colnames <- function(table, col_names) {
  test <- col_names %in% names(table)
  missing_cols <- col_names[!test]
  if (!all(test)) {
    stop("Input table lacks column(s) ",
         paste(paste0("'", missing_cols, "'"), collapse = ", "),
         ".")
  }
}

#' Check that the input table does not contain any missing values.
#' @param table A data.table passed to a scan statistic function.
#' @param col_names A character vector of column names; these columns in the
#'    table must not have any missing values.
#' @keywords internal
validate_values <- function(table, col_names) {
  if (any(is.na(table[, col_names, with = FALSE]))) {
    stop("The columns ",
         paste(paste0("'", col_names, "'"), collapse = ", "),
         " of the input table must not contain any missing values.")
  }
}

#' Check that the zones argument is a list of integer or factor vectors.
#' @param zones Should be a list of integer or factor vectors.
#' @keywords internal
validate_zones <- function(zones) {
  if (class(zones) != "list") {
    stop("The argument 'zones' must be a list of integer or factor vectors.")
  }
  if (!all(vapply(zones, class, character(1)) %in% c("integer", "factor"))) {
    stop("The argument 'zones' must be a list of integer or factor vectors.")
  }
}

#' Check that input to scanstatistic function is valid.
#' @param table A data.table passed to a scan statistic function.
#' @param zones Should be a list of integer or factor vectors.
#' @param col_names A character vector of column names; these columns in the
#'    table must not have any missing values.
#' @keywords internal
validate_scan <- function(table, zones, col_names) {
  validate_colnames(table, col_names)
  validate_values(table, col_names)
  validate_zones(zones)
}