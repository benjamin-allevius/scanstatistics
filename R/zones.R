# Functions in this file:
#   coords_to_knn
#   dist_to_knn
#   knn_zones
#   closest_subsets
#   flexible_zones
#   connected_neighbors
#   if_connected
#   is_connected
#   connected_to


#' Get the k nearest neighbors for each location, given its coordinates.
#' 
#' Get the k nearest neighbors for each location, including the location itself.
#' This function calls \code{\link[stats]{dist}}, so the options for the 
#' distance measure used is the same as for that one. Distances are calculated 
#' between rows.
#' @param k The number of nearest neighbors, counting the location itself.
#' @inheritParams stats::dist
#' @return An integer matrix of the \eqn{k} nearest neighbors for each location. 
#'    Each row corresponds to a location, with the first element of each row 
#'    being the location itself. Locations are encoded as integers.
#' @importFrom stats dist
#' @export
#' @examples 
#' x <- matrix(c(0, 0,
#'               1, 0,
#'               2, 1,
#'               0, 4,
#'               1, 3),
#'             ncol = 2, byrow = TRUE)
#' plot(x)
#' coords_to_knn(x)
coords_to_knn <- function(x, 
                          k = min(10, nrow(x)), 
                          method = "euclidean", 
                          p = 2) {
  dist_to_knn(dist(x, method = "euclidean", diag = T, upper = T, p = p), k = k)
}

#' Given a distance matrix, find the \eqn{k} nearest neighbors.
#' 
#' Given a distance matrix, calculate the \eqn{k} nearest neighbors of each 
#' location, including the location itself. The matrix should contain only zeros 
#' on the diagonal, and all other elements should be positive. 
#' @param x A (square) distance matrix. Elements should be non-negative and the 
#'    diagonal zeros, but this is not checked.
#' @inheritParams coords_to_knn
#' @return A matrix of integers, row \eqn{i} containing the \eqn{k} nearest 
#'    neighbors of location \eqn{i}, including itself.
#' @export
#' @examples 
#' x <- matrix(c(0, 0,
#'               1, 0,
#'               2, 1,
#'               0, 4,
#'               1, 3),
#'             ncol = 2, byrow = TRUE)
#' d <- dist(x, diag = TRUE, upper = TRUE)
#' dist_to_knn(d, k = 3)
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

#' Find the increasing subsets of \eqn{k} nearest neighbors for all locations.
#' 
#' Returns the set of increasing nearest neighbor sets for all locations, as
#' a list of integer vectors. That is, for each location the list returned 
#' contains one vector containing the location itself, another containing the 
#' location and its nearest neighbor, and so on, up to the vector containing the 
#' location and its \eqn{k-1} nearest neighbors. 
#' @param k_nearest An integer matrix of with \eqn{k} columns and as many rows
#'    as locations. The first element of each row is the integer encoding the
#'    location (and equal to the row number); the following elements are the 
#'    \eqn{k-1} nearest neighbors in ascending order of distance.
#' @importFrom plyr alply
#' @importFrom magrittr %>%
#' @return A list of integer vectors.
#' @export
#' @examples 
#' nn <- matrix(c(1L, 2L, 4L, 3L, 5L,
#'                2L, 1L, 3L, 4L, 5L, 
#'                3L, 2L, 4L, 1L, 5L,
#'                4L, 1L, 2L, 3L, 5L,
#'                5L, 3L, 4L, 2L, 1L),
#'                ncol = 5, byrow = TRUE)
#' knn_zones(nn[, 1:3])
knn_zones <- function(k_nearest) {
  k_nearest %>%
    unname %>%
    alply(.margins = 1, .fun = closest_subsets, .expand = F) %>%
    unlist(recursive = FALSE) %>%
    unname %>%
    unique
}


#' Set of increasing sets from left to right of input vector.
#' 
#' Returns a set (list) of the increasing sets (integer vectors) of the input 
#' vector \code{v}, in the sense that the first set contains the first element 
#' of \code{v}, the second set the first and second elements of \code{v}, and so 
#' on.
#' @param v An integer vector. Meant to represent the \eqn{k} nearest neighbors
#'    of a location, the first element being the integer identifier of the 
#'    location itself.
#' @return A list of the same length as the input. The first element of the list
#'    is v[1], the second is sort(v[1:2]), the third sort(v[1:3]), and so on.
#' @keywords internal
closest_subsets <- function(v) {
  lapply(seq_along(v), function(x) sort(v[seq(x)]))
}

#' Computes the flexibly shaped zones as in Tango (2005).
#' 
#' Given a matrix of \eqn{k} nearest neighbors and an adjacency matrix
#' for the locations involved, produces the set of flexibly shaped zones
#' as a list of integer vectors. The locations in these zones are all connected, 
#' in the sense that any location in the zone can be reached from another by 
#' traveling through adjacent locations within the zone.
#' @param k_nearest An integer matrix of the \eqn{k} nearest neighbors for each 
#'    location. Each row corresponds to a location, with the first element of 
#'    each row being the location itself. Locations should be encoded as 
#'    integers.
#' @inheritParams connected_to
#' @importFrom plyr alply
#' @importFrom magrittr %>%
#' @export
#' @return A list of integer vectors.
#' @references 
#'    Tango, T. & Takahashi, K. (2005), \emph{A flexibly shaped spatial scan 
#'    statistic for detecting clusters}, International Journal of Health 
#'    Geographics 4(1).
#' @examples 
#' A <- matrix(c(0,1,0,0,0,0,
#'               1,0,1,0,0,0,
#'               0,1,0,0,0,0,
#'               0,0,0,0,1,0,
#'               0,0,0,1,0,0,
#'               0,0,0,0,0,0), 
#'               nrow = 6, byrow = TRUE) == 1
#' nn <- matrix(as.integer(c(1,2,3,4,5,6,
#'                           2,1,3,4,5,6,
#'                           3,2,1,4,5,6,
#'                           4,5,1,6,3,2,
#'                           5,4,6,1,3,2,
#'                           6,5,4,1,3,2)),
#'                           nrow = 6, byrow = TRUE)
#' flexible_zones(nn, A)
flexible_zones <- function(k_nearest, adjacency_matrix) {
  k_nearest %>%
    unname %>%
    alply(.margins = 1, .fun = connected_neighbors, 
          adjacency_matrix = adjacency_matrix, .expand = F) %>%
    unlist(recursive = FALSE) %>%
    unname %>%
    unique %>%
    lapply(unlist)
}

#' Find the connected sets for a location and its \eqn{k} nearest neighbors.
#' 
#' Returns a \code{set} of \code{set}s, each set of the latter type containing
#' the location itself and zero or more of its neighbors, if they are connected.
#' @param neighbors A vector of neighbors to a location, the first element
#'    of the vector being the specific location, and the other elements its 
#'    other nearest neighbors. Locations should be encoded as integers.
#' @inheritParams connected_to
#' @importFrom sets set_power as.set set_union
#' @keywords internal
#' @return Returns a \code{set} of \code{set}s, each set of the latter type 
#'    containing the location itself and zero or more of its neighbors, if they 
#'    are connected.
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
