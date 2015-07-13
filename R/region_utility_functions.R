
#' Decide how many parts to partition the set of all regions into.
#' 
#' Decide how many parts to partition the set of all regions into, for applying
#' some function(s) over each part. Few parts can lead to memory problems for 
#' large data sets, while many parts lead to longer computation times.
#' @inheritParams partition_regions
auto_region_partition_size <- function(regions) {
  min(length(regions), floor(log(sum(vapply(regions, length, integer(1))))))
}

#' Partition a set of regions.
#' 
#' Partition a set of regions such that each part contains about the same
#' number of locations, when the number of locations in each region for the
#' part are summed over all regions in the part.
#' @param regions A \code{set} of regions, each region itself being a \code{set}
#'    containing locations.
#' @param n_parts An integer; the number of parts to split the \code{regions}
#'    into.
#' @return A list with two elements:
#'    \itemize{
#'      \item{partition} A list, each element of which is a \code{set} 
#'        containing one or more regions (\code{set} containing locations).
#'      \item{offsets} An integer vector containing offset numbers to the
#'        region numbering. For example, the first region in 
#'        \code{partition[i]} will have will be region number 
#'        \code{offset[i] + 1}.
#'    }
partition_regions <- function(regions, n_parts = min(10L, length(regions))) {
  if (n_parts < 1 || n_parts %% 1 != 0) {
    stop("n_parts has to be a positive integer.")
  }
  if (n_parts > length(regions)) {
    stop("Can't partition set of all regions into more parts than the total ",
         "number of regions in it.")
  }
  n_parts <- as.integer(n_parts)
  
  # Turn into list in order to be able to access ranges of elements by index
  regs <- as.list(regions)
  
  # Decide partition by looking at cumulative sum of number of locations
  n_locations <- vapply(regions, length, integer(1))
  cs <- cumsum(n_locations)
  total_n <- cs[length(cs)]
  
  # Get breakpoints for where to partition regions
  ranges <- integer(n_parts)
  for (r in seq(n_parts - 1)) {
    ranges[r] <- sum(cs <= floor(total_n / n_parts) * r)
  }
  ranges[n_parts] <- length(regions)
  
  # offsets keep track of the region numbers
  offsets <- c(0L, ranges[-n_parts])
  
  region_partition <- list()
  os <- 0L
  for (r in ranges) {
    region_partition <- c(region_partition, 
                          sets::set(sets::as.set(regs[(os + 1L):r])))
    os <- r
  }
  list(partition = region_partition, offsets = offsets)
}

#' Applies a function over all regions containing the locations supplied.
#' 
#' Applies the function \code{f} to the data table formed by expanding the
#' \code{location_table} according to the regions in \code{region_partition}.
#' @param location_table A \code{data.table} with key column \code{location},
#'    and others which may be used by the supplied function.
#' @param region_partition A list as outputted by 
#'    \code{\link{partition_regions}}. Has two elements:
#'    \describe{
#'      \item{\code{partition}}{A list, each element of which is a \code{set} 
#'        containing one or more regions (\code{set} containing locations).}
#'      \item{\code{offsets}}{An integer vector containing offset numbers to the
#'         region numbering. For example, the first region in 
#'         \code{partition[i]} will have will be region number 
#'         \code{offsets[i] + 1}.}
#'    }
#' @param f A function to apply after expanding \code{location_table} by the
#'    regions in a given element of \code{region_partition$partition}.
#' @param keys A character vector to set the key columns of the expanded 
#'    region-location table by before applying the function \code{f}.
#' @return A \code{data.table}, containing the results of applying the supplied
#'    function over all regions.
#' @importFrom foreach %do%
region_apply <- function(location_table, region_partition, f, keys = NULL) {
  foreach::foreach(regions_in_part = region_partition$partition, 
                   offset = region_partition$offsets,
                   .combine = rbind,
                   .packages = c("data.table", "magrittr")) %do% {
    locations <- unique(unlist(regions_in_part))
    region_table <- region_table_creator(regions_in_part, 
                                         keys = c("location"), 
                                         offset = offset)
    merge(location_table[location %in% locations, ],
          region_table,
          by = "location",
          allow.cartesian = TRUE) %>% {
            setkeyv(., keys)
            .
          } %>% f
  }
}


#' Creates a new \code{data.table} from a table containing locations,
#' adding a column for region.
#' 
#' Takes a \code{data.table} with containing column \code{location} and 
#' possibly other columns, and creates a new \code{data.table} with 
#' a column for region added to the columns in the supplied table, 
#' according to the regions in the supplied list of regions.
#' The key colums of the resulting \code{data.table} can be specified.
#' @param locations_etc A \code{data.table} with column \code{location}
#'    and other columns (but none for \code{region}).
#' @param regions A list of regions, elements being vectors of locations.
#' @param keys Character vector of one or more column names; these columns
#'    are set as key columns in the output \code{data.table}.
#' @return A new \code{data.table} with a column for \code{region} added
#'    to the supplied table of locations etc. (not modified).
#' @examples
#' \dontrun{
#' locs_etc <- table_creator(list(location = 1:2, duration = 1:3, stream = 1:2))
#' region_joiner(locs_etc, list(1, 2, 1:2), keys = c("duration", "region"))
#' }
region_joiner <- function(locations_etc, regions, keys = c("region")) {
  region_table_creator(regions, keys = "location")[
    locations_etc, allow.cartesian = TRUE][, .SD, keyby = keys]
}

#' Converts a list of regions to a \code{data.table} of regions and locations.
#' 
#' Supply a list of regions, with each element of the list being a vector
#' of locations. If the list is named, the output \code{data.table}
#' will have these names for the regions. Else, the regions are labeled
#' by integers from 1 to the length of the region list.
#' @param regions A list of regions, elements being vectors of locations.
#' @param keys Character vector of one or more column names which is passed 
#'    to \code{\link[data.table]{setkey}}.
#' @param offset An integer to offset the region numbering by, in case names are
#'    not used for the regions, and you want the region count to start at
#'    \code{offset} + 1.
#' @examples 
#' \dontrun{
#' region_table_creator(list(1L, 2L, 1:2))
#' region_table_creator(sets::set(sets::set(1L), 
#'                      sets::set(2L), sets::as.set(1:2)))
#' region_table_creator(list(1L, 2L, 1:2), keys = "location")
#' region_table_creator(list(1L, 2L, 1:2), keys = "region")
#' region_table_creator(list(a = "x", b = "y", c = c("x", "y")))
#' }
region_table_creator <- function(regions, keys = NULL, offset = 0L) {
  region_names <- names(regions)
  if (is.null(region_names)) {
    region_names <- seq_along(regions) + offset
  }
  data.table(location = unlist(regions, use.names = FALSE),
             region = rep(region_names, 
                          vapply(regions, length, integer(1))),
             key = keys)
}