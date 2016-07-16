# Functions in this file:
#   zone_sum
#   cumsum_duration


#' Sum columns over all location in each zone, for each duration.
#' 
#' @param table A \code{data.table} with columns \code{zone, location, duration} 
#'    and those given in the argument \code{sumcols}.
#' @param sumcols Character vector of column names, the columns to be summed 
#'    over each zone and duration.
#' @keywords internal
zone_sum <- function(table, sumcols) {
  e <- parse(text = paste0("list(", 
                           paste0(sumcols, " = sum(", sumcols, ")", 
                                  collapse = ", "), 
                           ")"))
  table[, eval(e), by = .(zone, duration)]
}


#' Calculates the cumulative sum of columns over duration, by other columns.
#' 
#' This function calculates the cumulative sum over the (column) duration for
#' one or more columns, by some other columns.
#' @param table A \code{data.table} with column \code{duration} and other 
#'    columns matching the other arguments to the function.
#' @param sumcols A character vector of the columns to calculate the cumulative
#'    sum for.
#' @param bycols A character vector of the columns to calculate the cumulative
#'    sum by (i.e. as the argument \code{by} to a \code{data.table}).
#' @return A \code{data.table} with the cumulative sums over duration.
#' @keywords internal
cumsum_duration <- function(table, sumcols, bycols) {
  e <- parse(text = paste0("list(", 
                           "duration = duration, ",
                           paste0(sumcols, " = cumsum(", sumcols, ")", 
                                  collapse = ", "), 
                           ")"))
  b <- parse(text = paste0("list(", 
                           paste0(bycols, collapse = ", "), 
                           ")"))
  table[, eval(e), by = eval(b)]
}