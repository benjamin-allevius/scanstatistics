
`%notin%` <- function(a, b) !(a %in% b)

# get package names
get_package_names <- function() {
  gsub("package:", "", search()[grep("package:", search())])
}

#' Is the relative error between two numbers is less than the given tolerance?
#' 
#' Given two consecutive numbers in a sequence, return \code{TRUE} if the
#' relative change is positive but less than the given tolerance.
#' @param current A scalar; the most recent value of the sequence.
#' @param previous A scalar; the second most recent value of the sequence, or a
#'    reference value.
#' @param tol The tolerance, a positive scalar near zero.
#' @keywords internal
has_converged <- function(current, previous, tol = 0.01) {
  rel_change <- (current - previous) / abs(previous)
  rel_change > 0 && rel_change < tol
}

#' Flip a matrix upside down
#' @param x A matrix
#' @return A matrix, \code{x} with rows reversed.
#' @keywords internal
flipud <- function(x) {
  x[rev(seq_len(nrow(x))), , drop = FALSE]
}

#' Convert a long data frame to a wide matrix.
#' 
#' Convert a long data frame to a wide matrix, with time along the row dimension
#' and locations along the column dimension. Values in the matrix could be e.g.
#' the observed counts or the population.
#' @param df A data frame with at least 3 columns.
#' @param time_col Integer or string that specifies the time column.
#' @param location_col Integer or string that specifies the location column.
#' @param value_col Integer or string that specifies the value column.
#' @return A matrix with time on rows and locations on columns.
#' @importFrom dplyr select
#' @importFrom tidyr spread_
#' @importFrom magrittr %>%
#' @export
df_to_matrix <- function(df, time_col = 1, location_col = 2, value_col = 3) {
  if (is.numeric(time_col)) {
    time_col <- names(df)[time_col]
  }
  if (is.numeric(location_col)) {
    location_col <- names(df)[location_col]
  }
  if (is.numeric(value_col)) {
    value_col <- names(df)[value_col]
  }
  x <- df %>%
    as.data.frame %>%
    select(time_col, location_col, value_col) %>%
    spread_(key_col = location_col, value_col = value_col)
  times <- x[, 1]
  x <- as.matrix(x[, -1])
  attributes(x)$dimnames[[1]] <- times
  return(x)
}


# Clean up when package is unloaded.
.onUnload <- function (libpath) {
  library.dynam.unload("scanstatistics", libpath)
}
