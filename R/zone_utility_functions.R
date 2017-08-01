# Functions in this file:
#   get_zone
#   powerset_zones

#' Extract a zone from the set of all zones.
#' 
#' Extract zone number \eqn{n} from the set of all zones.
#' @param n An integer; the number of the zone you wish to retrieve.
#' @param zones A list of integer vectors, representing the set of all zones.
#' @return An integer vector.
#' @export
#' @examples 
#' zones <- list(1L, 2L, 3L, 1:2, c(1L, 3L), c(2L, 3L))
#' get_zone(4, zones)
get_zone <- function(n, zones) {
  if (n > length(zones) || n < 1) {
    stop("Zone not found.")
  }
  zones[[n]]
}

#' Creates a set of all non-empty subsets of the integers from 1 to \eqn{n}.
#' 
#' Creates a list of all \eqn{2^(n-1)} non-empty subsets of the integers from 1 
#' to \eqn{n}.
#' @param n An integer larger than 0.
#' @return A list of integer vectors.
#' @importFrom utils combn
#' @keywords internal
powerset_zones <- function(n) {
  zones <- sets::set_power(seq_len(n)) - sets::set(sets::as.set(integer(0)))
  lapply(zones, FUN = function(x) unlist(as.list(x)))
}
