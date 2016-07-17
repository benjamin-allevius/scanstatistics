


#' Get the k nearest neighbors for each point, given its coordinates.
#' 
#' Get the k nearest neighbors for each point, including the point itself.
#' This function calls \code{\link[stats]{dist}},
#' so the options for the distance measure used is the same as for that one.
#' Distances are calculated between rows.
#' @param k The number of nearest neighbors, counting the point itself.
#' @inheritParams stats::dist
#' @return A matrix of integers, row \eqn{i} containing the \eqn{k} nearest 
#'    neighbors of point \eqn{i}, including itself.
#' @importFrom stats dist
#' @keywords internal
#' @export
coords_to_knn <- function(x, 
                          k = min(10, nrow(x)), 
                          method = "euclidean", 
                          p = 2) {
  dist_to_knn(dist(x, method = "euclidean", diag = T, upper = T, p = p), k = k)
}

#' Given a distance matrix, calculate \eqn{k} nearest neighbors.
#' 
#' Given a distance matrix, calculate the \eqn{k} nearest neighbors of each 
#' point, including the point itself. The matrix should contain only zeros on 
#' the diagonal, and all other elements should be positive. 
#' @param x A distance matrix. Elements should be non-negative and the diagonal
#'    zeros, but this is not checked.
#' @inheritParams coords_to_knn
#' @return A matrix of integers, row \eqn{i} containing the \eqn{k} nearest 
#'    neighbors of point \eqn{i}, including itself.
#' @keywords internal
dist_to_knn <- function(x, k = min(10, nrow(x))) {
  if (class(x) == "dist" && (!attr(x, "Diag") || !attr(x, "Upper"))) {
    stop("If x is a 'dist' object, it must have diag=TRUE and upper=TRUE")
  }
  t(apply(as.matrix(x), 2, order))[, seq(k)]
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
#' @param k_nearest A matrix or data frame with \eqn{k} columns, the \eqn{j}th 
#'    element of the \eqn{i}th row containg the \eqn{j}th nearest neighbor of 
#'    point \eqn{i}, with \eqn{j = 1,\ldots,k}. Each point is defined to be its 
#'    own nearest neighbor, so the first column of row \eqn{i} should be (the 
#'    identifier) of point \eqn{i} itself.
#' @inheritParams plyr::alply
#' @importFrom sets set_union
#' @importFrom plyr alply
knn_zones <- function(k_nearest, .parallel = FALSE, .paropts = NULL) {
  Reduce(set_union, 
         alply(k_nearest, 
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
#' @param v A vector. Meant to represent the \eqn{k} nearest neighbors
#'    of a point, the first element being (an identifier) of the point itself.
#' @importFrom sets as.set
#' @keywords internal
closest_subsets <- function(v) {
  as.set(lapply(lapply(seq_along(v), function(x) v[seq(x)]), as.set))
}

#' Computes the flexibly shaped zones as in Tango (2005).
#' 
#' Given a matrix of \eqn{k} nearest neighbors and an adjacency matrix
#' for the locations involved, produces the set of flexibly shaped zones
#' (sets of locations). The zones in these sets are all connected,
#' in the sense that any location in the zone can be reached from another
#' by traveling through adjacent locations within the zone.
#' @param k_nearest A matrix of the \eqn{k} nearest neighbors for each location.
#'    Each row corresponds to a location, with the first element of each row
#'    being the location itself. Locations should preferably be given as 
#'    integers.
#' @inheritParams connected_to
#' @inheritParams plyr::alply
#' @importFrom sets set_union as.set
#' @importFrom plyr alply
flexible_zones <- function(k_nearest, 
                             adjacency_matrix,
                             .parallel = FALSE, 
                             .paropts = NULL) {
  Reduce(set_union,
         as.set(alply(k_nearest,
                      .margins = 1,
                      .fun = connected_neighbors,
                      adjacency_matrix = adjacency_matrix,
                      .parallel = .parallel,
                      .paropts = .paropts)))
}

#' Returns the connected sets for a location and its \eqn{k} nearest neighbors,
#' including itself.
#' 
#' Returns a \code{set} of \code{set}s, each set of the latter type containing
#' the location itself and zero or more of its neighbors, if they are connected.
#' @param neighbors A vector of neighbors to a location, the first element
#'    of the vector being the specific location, and the other elements its 
#'    other nearest neighbors. Locations should preferably be integers.
#' @inheritParams connected_to
#' @importFrom sets set_power as.set set_union
#' @keywords internal
connected_neighbors <- function(neighbors, adjacency_matrix) {
  location <- neighbors[1]
  its_neighbors <- neighbors[-1]
  pset <- set_power(as.set(its_neighbors)) - sets::set(sets::set())
  set_union(sets::set(sets::set(location)),
            as.set(lapply(pset, 
                          if_connected, 
                          location = location,
                          adjacency_matrix = adjacency_matrix))) - 
    sets::set(sets::set())
}


#' Return a set of the location and its neighbors if they are connected,
#' else return the empty set.
#' 
#' If the location and its neighbors, not including itself, are connected,
#' then return the set containing the location and its neighbors;
#' otherwise, return the empty set
#' @param distinct_neighbors A \code{set} containing the neighboring locations
#'    to the given location, not including the location itself.
#' @param location A location, preferably given as an integer.
#' @inheritParams connected_to
#' @return A \code{set} of the given location and the neighbors if they are
#'    connected, else returns the empty set.
#' @importFrom sets set_union
#' @keywords internal
if_connected <- function(distinct_neighbors, location, adjacency_matrix) {
  if (is_connected(distinct_neighbors, location, adjacency_matrix)) {
    return(set_union(sets::set(location), distinct_neighbors))
  } else {
    return(sets::set())
  }
}


#' Returns TRUE if the neighboring locations are connected to the given 
#' location, FALSE if not.
#' 
#' @param neighbor_locations A \code{set} of neighboring locations to the given
#'    location; these neighbors do not include the given location itself.
#' @param location A location, preferably given as an integer.
#' @inheritParams connected_to
#' @return Boolean: is the neighbors connected to the given location?
#' @keywords internal
#' @importFrom sets set_is_empty
is_connected <- function(neighbor_locations, location, adjacency_matrix) {
  Z_0 <- sets::set(location)
  Z_1 <- neighbor_locations
  while (TRUE) {
    Z_0 <- connected_to(Z_0, Z_1, adjacency_matrix)
    if (set_is_empty(Z_0)) {
      return(FALSE)
    }
    Z_1 <- Z_1 - Z_0
    if (set_is_empty(Z_1)) {
      return(TRUE)
    }
  }
}

#' Return those elements in the second set which are connected to those in the
#' first.
#' 
#' Return those elements in the second set \eqn{Z_1} which are connected to 
#' those in the first set \eqn{Z_0}, according to the adjacency matrix.
#' @param Z_0 A set of locations, given as integers.
#' @param Z_1 A set of locations, given as integers.
#' @param adjacency_matrix A boolean matrix, with element \eqn{(i,j)} set 
#'    to TRUE if location \eqn{j} is adjacent to location \eqn{i}.
#' @return A set, possibly empty, containing those locations in \eqn{Z_1}
#'         that are connected to any of the locations in \eqn{Z_0}.
#' @importFrom sets set_union
#' @keywords internal
connected_to <- function(Z_0, Z_1, adjacency_matrix) {
  connected <- sets::set()
  for (loc in Z_1) {
    if (any(Z_0 %in% which(adjacency_matrix[loc, ]))) {
      connected <- set_union(connected, loc)
    }
  }
  return(connected)
}