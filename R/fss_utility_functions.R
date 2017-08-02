
#' Aggregate (sum) values over all data streams, and cumulatively over time.
#' @param A An array with three dimensions. Dimensions are: 
#'    \describe{
#'      \item{Dimension 1}{Time, ordered from most recent to most distant.}
#'      \item{Dimension 2}{Location, enumerated from 1 and up.}
#'      \item{Dimension 3}{Data stream, enumerated from 1 and up.}
#'    }  
#' @return A matrix with \code{dim(A)[1]} rows and \code{dim(A)[2]} columns.
#' @keywords internal
aggregate_per_location <- function(A) {
  apply(apply(A, 1:2, sum), 2, cumsum)
}

#' Aggregate (sum) values over all locations, and cumulatively over time.
#' @param A An array with three dimensions. Dimensions are: 
#'    \describe{
#'      \item{Dimension 1}{Time, ordered from most recent to most distant.}
#'      \item{Dimension 2}{Location, enumerated from 1 and up.}
#'      \item{Dimension 3}{Data stream, enumerated from 1 and up.}
#'    }  
#' @return A matrix with \code{dim(A)[1]} rows and \code{dim(A)[3]} columns.
#' @keywords internal
aggregate_per_stream <- function(A) {
  apply(apply(A, c(1, 3), sum), 2, cumsum)
}

#' Apply a function to each row of a matrix.
#' 
#' Apply a function to each row of a matrix. If the function returns a scalar,
#' then return a vector. If the function returns a vector, return a matrix with
#' the same number of rows.
#' @param A A matrix.
#' @param .f A function taking a vector as a first argument. This function 
#'    should preferably return a vector of the same length as the first 
#'    argument.
#' @param ... Other arguments passed to \code{.f}.
#' @return In case \code{.f} returns a vector of length \code{n}, a matrix with
#'    \code{nrow(A)} rows and \code{n} columns is returned. If \code{.f} returns
#'    a scalar, a vector with \code{nrow(A)} elements is returned.
#' @keywords internal
apply_rowwise <- function(A, .f, ...) {
  res <- apply(A, 1, .f, ...)
  if (is.null(dim(res))) {
    return(res)
  } else {
    return(t(res))
  }
}

#' Order row contents by priority for each timepoint.
#' 
#' Given a matrix of priority function values, return a matrix in which each 
#' row gives the column indices of the corresponding row in the priority matrix,
#' when that row has been sorted from highest to lowest priority value.
#' @param priority_matrix A numeric or integer matrix. Rows represent time (ordered from 
#'    most recent to most distant), columns represent e.g. locations or data
#'    streams (numbered from 1 and up). The element in row \eqn{i} and column
#'     \eqn{j} holds the priority of the \eqn{j}th location/data stream for 
#'     times \eqn{1,\ldots,i}.
#' @return A matrix of the same size as the input. On each row, column indices
#'    are given in order of priority.
#' @keywords internal
prioritize_cols <- function(priority_matrix) {
  
  # For each row (time), rank each value from smallest (rank 1) to largest.
  # When priority values are tied, columns indices with lower number go first.
  ranked_prios <- apply_rowwise(priority_matrix, 
                                function(x) order(x, rev(seq_along(x))))
  
  # For each row, replace the rank with the index of the corresponding column
  t(apply(ranked_prios, 1, function(x) rev(seq_len(ncol(ranked_prios))[x])))
}

#' Reorder rows by priority.
#' 
#' Reorder each row in the input matrix \code{A} by the column indices found
#' in the corresponding row of the matrix \code{priority_indices}.
#' @param A A matrix, containing e.g. counts or baselines. Rows represent time 
#'    (ordered from most recent to most distant), columns represent e.g. 
#'    locations or data streams (numbered from 1 and up).
#' @param priority_indices An integer matrix as output by 
#'    \code{\link{prioritize_cols}}.
#' @return An integer matrix of the same dimension as \code{A}.
#' @keywords internal
reorder_rows <- function(A, priority_indices) {
  t(sapply(seq_len(nrow(A)), function(t) A[t, priority_indices[t, ]]))
}

#' Order locations accorder to priority, then apply function.
#' @param .f A function taking a vector as first argument.
#' @param A A matrix, e.g. containing counts or baselines. Rows represent time 
#'    (ordered from most recent to most distant), columns represent locations 
#'    (numbered from 1 and up).
#' @param ... Further arguments passed to \code{.f}.
#' @return A matrix of the same dimension as \code{A}.
#' @keywords internal
prioritize_and_execute <- function(.f, A, prioritized_locations, ...) {
  apply_rowwise(reorder_rows(A, prioritized_locations), .f, ...)
}

#' Get either a location or stream subset.
#' @param lst A list.
#' @param d Integer; either 2 (locations) or 3 (streams).
#' @return A list.
#' @keywords internal
get_subset <- function(lst, subset, d = 3) {
  dims <- dim(lst[[1]])
  loc_subset <- seq_len(dims[2])
  stream_subset <- seq_len(dims[3])
  if (d == 2) {
    loc_subset <- subset
  } else if (d == 3) {
    stream_subset <- subset
  } else {
    stop("d must be either 2 or 3")
  }
  lapply(lst, function(x) x[ , loc_subset, stream_subset, drop = FALSE])
}


#' For each array in a list, sum over dimension \code{d}.
#' @param lst List containing 3-dimensional arrays. The second dimension should
#'    represent locations, and the third dimension should represent data 
#'    streams.
#' @param subset An integer vector holding a subset of indices.
#' @param d An integer; the array dimension to sum over.
#' @return A list containing matrices. This list has the same names as the input
#'    list.
#' @keywords internal
sum_over_subset <- function(lst, subset, d = 3) {
  lapply(get_subset(lst, subset, d), function(x) apply(x, (1:3)[-d], sum))
}


# Used for priority functions --------------------------------------------------

#' Cumulatively sum, order by priority, then sum again.
#' 
#' For a given matrix, form a new matrix by cumulatively summing over its 
#' columns. Reorder the row contents in this matrix by the corresponding indices 
#' in the given priority matrix. Then cumulatively sum over each row, yielding
#' a new matrix of the same size as the input matrix.
#' @param A A numeric or integer matrix.
#' @param priority_indices An integer matrix as output by 
#'    \code{\link{prioritize_cols}}.
#' @keywords internal
sum_reorder_sum <- function(A, priority_indices) {
  t(apply(reorder_rows(apply(A, 2, cumsum), priority_indices), 1, cumsum))
}
