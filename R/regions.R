
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


flexible_regions <- function(k_nearest, 
                             adjacency_matrix,
                             .parallel = FALSE, 
                             .paropts = NULL) {
  connected_to <- pryr::partial(connected_to_full,
                                adjacency_matrix = adjacency_matrix)
  sets::as.set(plyr::alply(k_nearest,
                           .margins = 1,
                           .fun = connected_neighbors,
                           .parallel = .parallel,
                           .paropts = .paropts))
}


connected_neighbors <- function(neighbors) {
  location <- neighbors[1]
  its_neighbors <- neighbors[-1]
  pset <- sets::set_power(sets::as.set(its_neighbors)) - sets::set(sets::set())
  sets::set_union(sets::set(sets::set(location)),
                  sets::as.set(lapply(pset, 
                                      if_connected, 
                                      location = location)))
}

# If the location and its neighbors, not including itself, are connected,
# then return the set of the location and its neighbors,
# else return the empty set
if_connected <- function(distinct_neighbors, location) {
  if (is_connected(distinct_neighbors, location)) {
    return(sets::set_union(sets::set(location),
                           distinct_neighbors))
  } else {
    return(sets::set())
  }
}

# is the set of the location and its neighbors connected?
is_connected <- function(neighbor_locations, location) {
  Z_0 <- sets::set(location)
  Z_1 <- neighbor_locations
  while (TRUE) {
    Z_0 <- connected_to(Z_0, Z_1)
    if (sets::set_is_empty(Z_0)) {
      return(FALSE)
    }
    Z_1 <- Z_1 - Z_0
    if (sets::set_is_empty(Z_1)) {
      return(TRUE)
    }
  }
}

# returns a set of the elements in Z_1 connected any of the elements in Z_0,
# according to the adjacency_matrix
# Element (i,j) of the adjacency_matrix is TRUE if i is adjacent to j
# with no elements being ajacent to themselves (so element (i,i) is FALSE)
# assumes locations are integers
connected_to_full <- function(Z_0, Z_1, adjacency_matrix) {
  connected <- sets::set()
  for (loc in Z_1) {
    if (any(Z_0 %in% which(adjacency_matrix[loc, ]))) {
      connected <- sets::set_union(connected, loc)
    }
  }
  return(connected)
}