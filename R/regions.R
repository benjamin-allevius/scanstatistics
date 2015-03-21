
#' Get the k nearest neighbors for each point.
#' 
#' Get the k nearest neighbors for each point, including the point itself.
#' This function calls \code{\link[stats]{dist}},
#' so the options for the distance measure used is the same as for that one.
#' 
#' @param k The number of nearest neighbors, counting the point itself.
#' @inheritParams stats::dist
k_nearest_neighbors <- function(x, 
                                k = min(10, nrow(x)), 
                                method = "euclidean", 
                                p = 2) {
  t(apply(as.matrix(dist(x, method = "euclidean", diag = T, upper = T, p = 2)),
          2, order))[, seq(k)]
}


# k_nearest is a matrix where each row corresponds to a location,
# with element j of row i corresponding to the jth closest location
# to location i. Each location is its own 1st closest location,
# so element i of the first column is location i itself.

#' Increasing subsets of \eqn{k} nearest neighbors.
#' 
#' Returns the set of increasing nearest neighbor sets for all points.
#' That is, for each point the set returned contains a set containing
#' the point itself, a set containing the point and its nearest neighbor,
#' and so on, up to the set containing the point and its \eqn{k-1} nearest 
#' neighbors. The set returned contains no duplicates.
#' 
#' @param k_nearest A matrix or data frame with \eqn{k} columns,
#'        the \eqn{j}th element of the \eqn{i}th row containg the \eqn{j}th
#'        nearest neighbor of point \eqn{i}, with \eqn{j = 1,\ldots,k}.
#'        Each point is defined to be its own nearest neighbor,
#'        so the first column of row \eqn{i} should be (the identifier)
#'        of point \eqn{i} itself.
#' @inheritParams plyr::alply
regions_upto_k <- function(k_nearest, .parallel = FALSE, .paropts = NULL) {
  Reduce(sets::set_union, 
         plyr::alply(k_nearest, 
                     .margins = 1, 
                     .fun = closest_subsets, 
                     .expand = F, 
                     .parallel = FALSE, 
                     .paropts = NULL))
}


#' Set of increasing sets from left to right of input vector.
#' 
#' Returns a set of the increasing sets of the input vector \code{v},
#' in the sense that the first set contains the first element of \code{v},
#' the second set the first and second elements of \code{v}, and so on.
#' 
#' @param v A vector. Meant to represent the \eqn{k} nearest neighbors
#'        of a point, the first element being (an identifier) of the point
#'        itself.
closest_subsets <- function(v) {
  sets::as.set(lapply(lapply(seq_along(v), function(x) v[seq(x)]), 
                      sets::as.set))
}
