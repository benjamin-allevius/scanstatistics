

# log-sum-exp trick: avoids arithmetic underflow and overflow
logsumexp <- function(x) {
    A <- max(x)
    A + log(sum(exp(x - A)))
}

#' Create a \code{data.table} with all combinations of the supplied variables.
#' 
#' @param col_list A named list of the variables you want 
#'        in your \code{data.table}.
#' @return A \code{data.table} with all combinations of the variables
#'         supplied in \code{col_list}.
#' @examples
#' cols <- list(location = 1:2, 
#'              time = 0:2, 
#'              stream = 1:2)
#' table_creator(cols)
table_creator <- function(col_list) {
  data.table(do.call(expand.grid, col_list))
}


#' Converts a list of regions to a \code{data.table} of regions and locations.
#' 
#' Supply a list of regions, with each element of the list being a vector
#' of locations. If the list is named, the output \code{data.table}
#' will have these names for the regions. Else, the regions are labeled
#' by integers from 1 to the length of the region list.
#' 
#' @param regions A list of regions, elements being vectors of locations.
#' @examples 
#' region_table_creator(list(1, 2, 1:2))
region_table_creator <- function(regions) {
  
  region_names <- names(regions)
  if (is.null(region_names)) {
    region_names <- seq_along(regions)
  }
  
  data.table(location = unlist(regions, use.names = FALSE),
             region = rep(region_names, 
                          vapply(regions, length, integer(1))),
             key = "location")
}

#' Add a column for \code{region} to a \code{data.table} containing 
#' \code{location}s.
#' 
#' Takes a \code{data.table} with key column \code{location} and preferably
#' other columns, and creates a new \code{data.table} with a column for region 
#' added to the columns in the supplied table, according to the regions
#' in the supplied list of regions.
#' @param locations_etc A \code{data.table} with column \code{location}
#'        and other columns (but none for \code{region}).
#' @param regions A list of regions, elements being vectors of locations.
#' @return A new \code{data.table} with a column for \code{region} added
#'         to the supplied table of locations etc. (not modified).
#' @examples
#' locs_etc <- table_creator(list(location = 1:2, time = 0:2, stream = 1:2))
#' region_joiner(locs_etc, list(1, 2, 1:2))
region_joiner <- function(locations_etc, regions) {
  merge(locations_etc, region_table_creator(regions),
        by = "location", allow.cartesian = TRUE)
}
