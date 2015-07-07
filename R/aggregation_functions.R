
#' Sum columns over all locations in per region and duration.
#' 
#' 
region_sum <- function(table, sumcols) {
  e <- parse(text = paste0("list(", 
                           paste0(sumcols, " = sum(", sumcols, ")", 
                                  collapse = ", "), 
                           ")"))
  table[, eval(e), by = .(region, duration)]
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