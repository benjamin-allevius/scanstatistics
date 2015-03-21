

# get package names
get_package_names <- function() {
  gsub("package:", "", search()[grep("package:", search())])
}

# May over/underflow
logsumexp <- function(x) {
  log(sum(exp(x)))
}

# log-sum-exp trick: avoids arithmetic underflow and overflow
logsumexp_stable <- function(x) {
    A <- max(x)
    A + log(sum(exp(x - A)))
}

haskeys <- function(data_table, keys) {
  keys %in% getkeys(data_table)
}

# Get the keys from a data.table
getkeys <- function(data_table) {
  attributes(data_table)$sorted
}

# Check if the first few keys of the supplied table
# match those in keys
first_keys_are_equal <- function(data_table, keys) {
  table_keys <- getkeys(data_table)
  if (is.null(table_keys) || length(keys) > length(table_keys)) {
    return(FALSE)
  } else {
    return(any(keys == table_keys[1:length(table_keys)]))
  }
}


#' Create a \code{data.table} with all combinations of the supplied variables.
#' 
#' @param col_list A named list of the variables you want 
#'        in your \code{data.table}.
#' @param key Character vector of one or more column names which is passed 
#'        to code{\link[data.table]{setkey}}.
#' @return A \code{data.table} with all combinations of the variables
#'         supplied in \code{col_list}.
#' @examples
#' cols <- list(location = 1:2, time = 0:2, stream = 1:2)
#' table_creator(cols)
table_creator <- function(col_list, key = NULL) {
  data.table(do.call(expand.grid, col_list),
             key = key)
}


#' Converts a list of regions to a \code{data.table} of regions and locations.
#' 
#' Supply a list of regions, with each element of the list being a vector
#' of locations. If the list is named, the output \code{data.table}
#' will have these names for the regions. Else, the regions are labeled
#' by integers from 1 to the length of the region list.
#' 
#' @param regions A list of regions, elements being vectors of locations.
#' @param key Character vector of one or more column names which is passed 
#'        to code{\link[data.table]{setkey}}.
#' @examples 
#' region_table_creator(list(1, 2, 1:2))
#' region_table_creator(sets::set(sets::set(1L), sets::set(2L), 
#'                                sets::as.set(1:2)))
#' region_table_creator(list(1, 2, 1:2), key = "location")
#' region_table_creator(list(1, 2, 1:2), key = "region")
#' region_table_creator(list(a = "x", b = "y", c = c("x", "y")))
region_table_creator <- function(regions, key = NULL) {
  
  region_names <- names(regions)
  if (is.null(region_names)) {
    region_names <- seq_along(regions)
  }
  
  data.table(location = unlist(regions, use.names = FALSE),
             region = rep(region_names, 
                          vapply(regions, length, integer(1))),
             key = key)
}

#' Add a column for \code{region} to a \code{data.table} containing 
#' \code{location}s.
#' 
#' Takes a \code{data.table} with containing column \code{location} and 
#' preferably other columns, and creates a new \code{data.table} with 
#' a column for region added to the columns in the supplied table, 
#' according to the regions in the supplied list of regions.
#' 
#' @param locations_etc A \code{data.table} with column \code{location}
#'        and other columns (but none for \code{region}).
#' @param regions A list of regions, elements being vectors of locations.
#' @param key Character vector of one or more column names which is passed 
#'        to code{\link[data.table]{setkey}}.
#' @return A new \code{data.table} with a column for \code{region} added
#'         to the supplied table of locations etc. (not modified).
#' @examples
#' locs_etc <- table_creator(list(location = 1:2, time = 0:2, stream = 1:2))
#' region_joiner(locs_etc, list(1, 2, 1:2))
region_joiner <- function(locations_etc, regions) {
  merge(x = region_table_creator(regions, key = "location"), 
        y = locations_etc,
        by = "location", allow.cartesian = TRUE)
}


#' Applies a function over all regions containing the locations supplied.
#' 
#' Applies the function \code{f} to the data table formed by expanding the
#' \code{location_table} according to the regions in \code{region_partition}.
#'
#' @param f The function to apply.
#' @param region_partition A \code{list} of \code{list}s. 
#'        The elements of each sublist are regions, 
#'        which are \code{vector}s of locations contained in the region.
#' @param location_table A \code{data.table} with key column \code{location},
#'        and other which may be used by the supplied function \code{f}.
region_apply <- function(f, region_partition, location_table) {
  foreach(regions_in_part = region_partition, 
          .combine = rbind,
          #           .packages = c("data.table", "magrittr"),
          .export = ls(envir = globalenv())) %dopar% {
            
            locations <- unique(unlist(regions))
            
            region_table <- region_table_creator(regions_in_part, 
                                                 key = c("location"))
            
            merge(location_table[location %in% locations, ],
                  region_table,
                  allow.cartesian = TRUE) %>% f
          }
}