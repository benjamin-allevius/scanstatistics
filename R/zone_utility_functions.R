# Functions in this file:
#   partition_zones
#   zone_apply
#   join_zones
#   zone_table_creator
#   get_zone
#   powerset_zones

#' Partition a set of zones.
#' 
#' Partition a set of zones such that each part contains about the same number 
#' of locations, when the number of locations in each zone for the part are 
#' summed over all zones in the part.
#' @param zones A \code{set} of zones, each zone itself being a \code{set}
#'    containing locations. Locations should be encoded as integers.
#' @param n_parts An integer; the number of parts to split the \code{zones} 
#'    into.
#' @return A list with two elements:
#'    \itemize{
#'      \item{partition} A list, each element of which is a \code{set} 
#'        containing one or more zones (\code{set} containing locations).
#'      \item{offsets} An integer vector containing offset numbers to the
#'        zone numbering. For example, the first zone in 
#'        \code{partition[i]} will have will be zone number 
#'        \code{offset[i] + 1}.
#'    }
#' @importFrom sets as.set
#' @keywords internal
partition_zones <- function(zones, n_parts = min(10L, length(zones))) {
  if (n_parts < 1 || n_parts %% 1 != 0) {
    stop("n_parts has to be a positive integer.")
  }
  if (n_parts > length(zones)) {
    stop("Can't partition set of all zones into more parts than the total ",
         "number of zones in it.")
  }
  n_parts <- as.integer(n_parts)
  
  # Turn into list in order to be able to access ranges of elements by index
  regs <- as.list(zones)
  
  # Decide partition by looking at cumulative sum of number of locations
  n_locations <- vapply(zones, length, integer(1))
  cs <- cumsum(n_locations)
  total_n <- cs[length(cs)]
  
  # Get breakpoints for where to partition zones
  ranges <- integer(n_parts)
  for (r in seq(n_parts - 1)) {
    ranges[r] <- sum(cs <= floor(total_n / n_parts) * r)
  }
  ranges[n_parts] <- length(zones)
  
  # offsets keep track of the zone numbers
  offsets <- c(0L, ranges[-n_parts])
  
  zone_partition <- list()
  os <- 0L
  for (r in ranges) {
    zone_partition <- c(zone_partition, sets::set(as.set(regs[(os + 1L):r])))
    os <- r
  }
  list(partition = zone_partition, offsets = offsets)
}

# #' Applies a function over all zones containing the locations supplied.
# #' 
# #' Applies the function \code{f} to the data table formed by expanding the
# #' \code{location_table} according to the zones in \code{zone_partition}.
# #' @param location_table A \code{data.table} with key column \code{location},
# #'    and others which may be used by the supplied function.
# #' @param zone_partition A list as outputted by \code{\link{partition_zones}}. 
# #'    Has two elements:
# #'    \describe{
# #'      \item{\code{partition}}{A list, each element of which is a \code{set} 
# #'        containing one or more zones (\code{set} containing locations).}
# #'      \item{\code{offsets}}{An integer vector containing offset numbers to the
# #'         zone numbering. For example, the first zone in \code{partition[i]} 
# #'         will have will be zone number \code{offsets[i] + 1}.}
# #'    }
# #' @param f A function to apply after expanding \code{location_table} by the
# #'    zones in a given element of \code{zone_partition$partition}.
# #' @param keys A character vector to set the key columns of the expanded 
# #'    zone-location table by before applying the function \code{f}.
# #' @return A \code{data.table}, containing the results of applying the supplied
# #'    function over all zones.
# #' @importFrom foreach foreach %do%
# #' @keywords internal
# zone_apply <- function(location_table, zone_partition, f, keys = NULL) {
#   
#   foreach(zones_in_part = zone_partition$partition, 
#           offset = zone_partition$offsets,
#           .combine = rbind,
#           .packages = c("data.table", "magrittr")) %do% {
#     locations <- unique(unlist(zones_in_part))
#     zone_table <- zone_table_creator(zones_in_part, 
#                                      keys = c("location"), 
#                                      offset = offset)
#     merge(location_table[location %in% locations, ],
#           zone_table,
#           by = "location",
#           allow.cartesian = TRUE) %>% {
#             setkeyv(., keys)
#             .
#           } %>% f
#   }
# }


#' Creates a new \code{data.table} from a table containing locations,
#' adding a column for zone.
#' 
#' Takes a \code{data.table} with containing column \code{location} and 
#' possibly other columns, and creates a new \code{data.table} with 
#' a column for zone added to the columns in the supplied table, 
#' according to the zones in the supplied list of zones.
#' The key colums of the resulting \code{data.table} can be specified.
#' @param locations_etc A \code{data.table} with column \code{location}
#'    and other columns (but none for \code{zone}).
#' @param zones A list of zones, elements being vectors of locations.
#' @param keys Character vector of one or more column names; these columns
#'    are set as key columns in the output \code{data.table}.
#' @return A new \code{data.table} with a column for \code{zone} added
#'    to the supplied table of locations etc. (not modified).
#' @examples
#' \dontrun{
#' locs_etc <- table_creator(list(location = 1:2, duration = 1:3, stream = 1:2))
#' join_zones(locs_etc, list(1, 2, 1:2), keys = c("duration", "zone"))
#' }
#' @keywords internal
join_zones <- function(locations_etc, zones, keys = c("zone")) {
  zone_table_creator(zones, keys = "location")[
    locations_etc, allow.cartesian = TRUE][, .SD, keyby = keys]
}

#' Converts a list of zones to a \code{data.table} of zones and locations.
#' 
#' Supply a list of zones, with each element of the list being a vector
#' of locations. If the list is named, the output \code{data.table}
#' will have these names for the zones. Else, the zones are labeled
#' by integers from 1 to the length of the zone list.
#' @param zones A list of zones, elements being vectors of locations.
#' @param keys Character vector of one or more column names which is passed 
#'    to \code{\link[data.table]{setkey}}.
#' @param offset An integer to offset the zone numbering by, in case names are
#'    not used for the zones, and you want the zone count to start at
#'    \code{offset} + 1.
#' @examples 
#' \dontrun{
#' zone_table_creator(list(1L, 2L, 1:2))
#' zone_table_creator(sets::set(sets::set(1L), 
#'                      sets::set(2L), sets::as.set(1:2)))
#' zone_table_creator(list(1L, 2L, 1:2), keys = "location")
#' zone_table_creator(list(1L, 2L, 1:2), keys = "zone")
#' zone_table_creator(list(a = "x", b = "y", c = c("x", "y")))
#' }
#' @keywords internal
zone_table_creator <- function(zones, keys = NULL, offset = 0L) {
  zone_names <- names(zones)
  if (is.null(zone_names)) {
    zone_names <- seq_along(zones) + offset
  }
  data.table(location = unlist(zones, use.names = FALSE),
             zone = rep(zone_names, 
                          vapply(zones, length, integer(1))),
             key = keys)
}


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
  zones <- sets::set()
  for (k in 1:n) {
    ss <- combn(1:n, k)
    for (j in 1:ncol(ss)) {
      zones <- sets::set_union(zones, sets::set(sets::as.set(ss[, j])))
    }
  }
  lapply(zones, FUN = function(x) unlist(as.list(x)))
}