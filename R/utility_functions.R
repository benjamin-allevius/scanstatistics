
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

#' Convert a matrix to a data frame.
#' 
#' Convert a matrix to a data frame with columns time, location, and one more
#' containing the elements of the input matrix.
#' @param mat A matrix.
#' @param name The name of the third column in the output matrix.
#' @param locations If not \code{NULL}, a vector with the names of the 
#'    locations.
#' @param times If not \code{NULL}, a vector with the time points. If 
#'    \code{NULL}, the matrix is assumed to be ordered with time point 1 in the 
#'    first row.
#' @return A matrix with columns \code{time, location, name}, where \code{name}
#'    is specified in the input.
#' @keywords internal
matrix_to_df <- function(mat, name, locations = NULL, times = NULL) {
  if (!is.null(locations) && length(locations) != ncol(mat)) {
    stop("The number of locations must be equal to col(mat)")
  }
  if (!is.null(times) && length(times) != nrow(mat)) {
    stop("The number of times must be equal to nrow(mat)")
  }
  if (is.null(locations)) {
    locations <- seq_len(ncol(mat))
  }
  if (is.null(times)) {
    times <- seq_len(nrow(mat))
  }
  
  df <- data.frame(time = rep(times, ncol(mat)),
                   location = rep(locations, each = nrow(mat)),
                   V1 = as.vector(mat))
  names(df)[3] <- name
  return(df)
}

#' Run a scan statistic analysis.
#' 
#' Run a scan statistic analysis with the given scan statistic and arguments.
#' @param scanstat A scan statistic function.
#' @param args A named list of arguments to be passed to \code{scanstat}.
#' @return A list with components
#'    \describe{
#'      \item{observed}{The table of observed statistics.}
#'      \item{simulated}{The table of simulated statistics.}
#'      \item{MC_pvalue}{The Monte Carlo P-value of the scan statistic.}
#'      \item{Gumbel_pvalue}{The Gumbel P-value of the scan statistic.}
#'    }
#' @keywords internal
run_scan <- function(scanstat, args) {
  scan <- do.call(scanstat, args)
  
  # Extract the most likely cluster (MLC)
  scan$observed %<>% arrange(-score)
  MLC <- scan$observed[1, ]
  
  # Get P-values
  gumbel_pvalue <- NULL
  MC_pvalue <- NULL
  if (nrow(scan$simulated) > 0) {
    gumbel_pvalue <- gumbel_pvalue(MLC$score, scan$simulated$score, 
                                   method = "ML")$pvalue
    MC_pvalue <- mc_pvalue(MLC$score, scan$simulated$score)
  }
  
  return(list(observed = scan$observed,
              replicates = scan$simulated,
              MC_pvalue = MC_pvalue,
              Gumbel_pvalue = gumbel_pvalue))
}


# Clean up when package is unloaded.
.onUnload <- function (libpath) {
  library.dynam.unload("scanstatistics", libpath)
}
