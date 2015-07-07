

#' Extract value of scan statistic from per-window statistics.
#' 
#' This function extracts the value of the scan statistic, which is the maximum
#' of the statistics calculated for each spatial or space-time window.
#' @param table A \code{data.table} with a column \code{statistic}, which should
#'    correspond to the statistic calculated for each spatial or space-time 
#'    window (given as other columns).
#' @return The maximum value of the column \code{statistic}.
extract_scanstatistic <- function(table) {
  table[, max(statistic)]
}

#' Extract the most likely cluster (MCL) and the value of the scans statistic.
#' 
#' This function extracts the most likely cluster, which is the spatial or 
#' spatiotemporal window that corresponds to the scan statistic. It also returns
#' the value of the scan statistc.
#' @inheritParams extract_scanstatistic
#' @return The row of the input table with the highest value of the column 
#'    \code{statistic}.
extract_mcl <- function(table) {
  table[which.max(statistic), ]
}