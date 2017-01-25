
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

#' Order locations by priority for each timepoint.
#' 
#' Given a matrix of priority function values, order the locations by priority
#' for each time point.
#' @param priority_mat A numeric matrix. Rows represent time (ordered from most
#'    recent to most distant), columns represent locations (numbered from 1 and 
#'    up).
#' @return A matrix of the same size as the input. On each row, the location
#'    numbers are ordered by priority.
#' @keywords internal
prioritize_locations <- function(priority_mat) {
  
  # For each row (time), rank each value from smallest (rank 1) to largest
  ranked_prios <- apply_rowwise(priority_mat, order)
  
  # For each row, replace the rank with the number of the corresponding location
  t(apply(ranked_prios, 1, function(x) rev(seq_len(ncol(ranked_prios))[x])))
}

#' Reorder locations (rows) by priority.
#' @param A matrix, e.g. containing counts or baselines. Rows represent time 
#'    (ordered from most recent to most distant), columns represent locations 
#'    (numbered from 1 and up).
#' @param priod_locations An integer matrix.
#' @return An integer matrix of the same dimension as \code{A}.
#' @keywords internal
reorder_locations <- function(A, priod_locations) {
  t(sapply(seq_len(nrow(A)), function(t) A[t, priod_locations[t, ]]))
}

#' Order locations accorder to priority, then apply function.
#' @param A A matrix, e.g. containing counts or baselines. Rows represent time 
#'    (ordered from most recent to most distant), columns represent locations 
#'    (numbered from 1 and up).
#' @return A matrix of the same dimension as \code{A}.
#' @keywords internal
prioritize_and_execute <- function(.f, A, prioritized_locations) {
  apply_rowwise(reorder_locations(A, prioritized_locations), .f)
}


#' @param counts A matrix of counts, possibly aggregated. Rows indicate time,
#'    ordered from most recent to most distant. Columns indicate locations,
#'    enumerated from 1 and up.
#' @param baselines A matrix of expected counts, possibly aggregated. Dimensions
#'    are as for \code{counts}.
#' @param score_fun A function taking matrix or vector arguments \code{counts},
#'    \code{baselines}, and possibly others (given in \code{...}), all of the
#'    same dimension, and returning a matrix or vector of that dimension.
#' @param prio_fun A function taking matrix or vector arguments \code{counts},
#'    \code{baselines}, and possibly others (given in \code{...}), all of the
#'    same dimension, and returning a matrix or vector of that dimension.
#' @param ... Arguments passed to \code{score_fun} and/or \code{prio_fun}. 
#'    Must be of the same dimension as \code{counts} and \code{baselines}.
#' @return A list containing three elements:
#'    \describe{
#'      \item{score}{The highest score of all clusters.}
#'      \item{duration}{The duration of the score-maximizing cluster.}
#'      \item{locations}{An integer vector of the locations in the 
#'                       score-maximizing cluster.}
#' }
#' @keywords internal
algo1 <- function(counts, 
                  baselines, 
                  score_fun,
                  prio_fun,
                  ...) {
  
  # Compute location priorities and sort them thereafter
  prios <- prio_fun(counts, baselines, ...)
  priod_locations <- prioritize_locations(prios)
  
  # Reorder aggregates by priority and take cumulative sums
  counts <- apply_rowwise(reorder_locations(counts, 
                                            priod_locations), 
                          cumsum)
  baselines <- apply_rowwise(reorder_locations(baselines, 
                                               priod_locations), 
                             cumsum)
  
  # Reorder by priority and take cumulative sums
  args <- lapply(list(x = counts, b = baselines, ...),
                 prioritize_and_execute, 
                 .f = cumsum,
                 prioritized_locations = priod_locations)
                 
  scores <- do.call(score_fun, args)
  
#   # # Score each prioritized subset of locations
#   # scores <- matrix(NA, nrow(prios), ncol(prios))
#   # for (t in seq_len(nrow(prios))) {
#   #   for (i in seq_len(ncol(prios))) {
#   #     scores[t, i] <- score_fun(agg_c[t, i], agg_b[t, i])
#   #   }
#   # }
  max_score <- max(scores)
  maxer <- which(scores == max_score, arr.ind = TRUE)
  duration <- unname(maxer[1, 1])
  locations <- priod_locations[duration, seq_len(maxer[1, 2])]
  list(duration = duration,
       locations = locations,
       score = max_score)
}

FN_SA <- function(counts, 
                  baselines, 
                  score_fun = function(c, b) ifelse(c > b, c * log(c/b) + b - c, 0)) {
  stream_subsets <- lapply(sets::set_power(seq_len(dim(counts)[3])), unlist)[-1]
  res <- lapply(stream_subsets, 
                function(x) algo1(counts[ , , x], baselines[ , , x], score_fun))
  scores <- vapply(res, function(lst) lst[["score"]], numeric(1))
  maximizer <- which.max(scores)
  top_scoring <- res[[maximizer]]
  top_scoring$streams <- stream_subsets[[maximizer]]
  top_scoring
}
