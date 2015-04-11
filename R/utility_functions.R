

# get package names
get_package_names <- function() {
  gsub("package:", "", search()[grep("package:", search())])
}

# May over/underflow
logsumexp_unstable <- function(x) {
  log(sum(exp(x)))
}

# log-sum-exp trick: avoids arithmetic underflow and overflow
logsumexp <- function(x) {
    A <- max(x)
    A + log(sum(exp(x - A)))
}

haskeys <- function(data_table, keys) {
  keys %in% getkeys(data_table)
}

#' Get the keys from a data.table.
#' 
#' @param data_table A \code{data.table}.
#' @return NULL if the supplied \code{data.table} has no keys,
#'         else a character vector containing the keys.
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
#'        to \code{\link[data.table]{setkey}}.
#' @return A \code{data.table} with all combinations of the variables
#'         supplied in \code{col_list}.
#' @examples
#' \dontrun{
#' cols <- list(location = 1:2, time = 0:2, stream = 1:2)
#' table_creator(cols)
#' }
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
#'        to \code{\link[data.table]{setkey}}.
#' @param offset An integer to offset the region numbering by, in case names are
#'        not used for the regions, and you want the region count to start at
#'        \code{offset} + 1.
#' @examples 
#' \dontrun{
#' region_table_creator(list(1L, 2L, 1:2))
#' region_table_creator(sets::set(sets::set(1L), sets::set(2L), sets::as.set(1:2)))
#' region_table_creator(list(1L, 2L, 1:2), key = "location")
#' region_table_creator(list(1L, 2L, 1:2), key = "region")
#' region_table_creator(list(a = "x", b = "y", c = c("x", "y")))
#' }
region_table_creator <- function(regions, key = NULL, offset = 0L) {
  region_names <- names(regions)
  if (is.null(region_names)) {
    region_names <- seq_along(regions) + offset
  }
  data.table(location = unlist(regions, use.names = FALSE),
             region = rep(region_names, 
                          vapply(regions, length, integer(1))),
             key = key)
}

#' Creates a new \code{data.table} from a table containing locations,
#' adding a column for region.
#' 
#' Takes a \code{data.table} with containing column \code{location} and 
#' possibly other columns, and creates a new \code{data.table} with 
#' a column for region added to the columns in the supplied table, 
#' according to the regions in the supplied list of regions.
#' The key colums of the resulting \code{data.table} can be specified.
#' 
#' @param locations_etc A \code{data.table} with column \code{location}
#'        and other columns (but none for \code{region}).
#' @param regions A list of regions, elements being vectors of locations.
#' @param keys Character vector of one or more column names; these columns
#'        are set as key columns in the output \code{data.table}.
#' @return A new \code{data.table} with a column for \code{region} added
#'         to the supplied table of locations etc. (not modified).
#' @examples
#' \dontrun{
#' locs_etc <- table_creator(list(location = 1:2, time = 0:2, stream = 1:2))
#' region_joiner(locs_etc, list(1, 2, 1:2), keys = c("time", "region"))
#' }
region_joiner <- function(locations_etc, regions, keys = c("region")) {
  region_table_creator(regions, key = "location")[
    locations_etc, allow.cartesian = TRUE][, .SD, keyby = keys]
}

#' Applies a function over all regions containing the locations supplied.
#' 
#' Applies the function \code{f} to the data table formed by expanding the
#' \code{location_table} according to the regions in \code{region_partition}.
#' @param location_table A \code{data.table} with key column \code{location},
#'    and others which may be used by the supplied function.
#' @param region_partition A list as outputted by 
#'    \code{\link{partition_regions}}. Has two elements:
#'    \itemize{
#'      \item{partition} A list, each element of which is a \code{set} 
#'        containing one or more regions (\code{set} containing locations).
#'      \item{offsets} An integer vector containing offset numbers to the
#'         region numbering. For example, the first region in 
#'         \code{partition[i]} will have will be region number 
#'         \code{offset[i] + 1}.
#'    }
#' @param f A function to apply after expanding \code{location_table} by the
#'    regions in a given element of \code{region_partition$partition}.
#' @param key A character vector to set the key columns of the expanded 
#'    region-location table by before applying the function \code{f}.
#' @return A \code{data.table}, containing the results of applying the supplied
#'    function over all regions.
region_apply <- function(location_table, region_partition, f, key = NULL) {
  foreach(regions_in_part = region_partition$partition, 
          offset = region_partition$offsets,
          .combine = rbind,
          .packages = c("data.table", "magrittr"),
          .export = ls(as.environment("package:scanstatistics"))) %dopar% {
            
            locations <- unique(unlist(regions_in_part))
            region_table <- region_table_creator(regions_in_part, 
                                                 key = c("location"), 
                                                 offset = offset)
            merge(location_table[location %in% locations, ],
                  region_table,
                  by = "location",
                  allow.cartesian = TRUE) %>% {
                    setkeyv(., key)
                    .
                  } %>% f
          }
}

#' Add a column \code{duration} to a \code{data.table} with column \code{time}.
#' 
#' Adds a column \code{duration} to a \code{data.table} with column \code{time},
#' the most recent time given a duration of 1, the second most recent time 
#' given a duration of 2, and so on. This function \strong{modifies} the input
#' table.
#' 
#' @param d A \code{data.table} containing at least a column \code{time},
#'        which is sortable. For example, could be POSIXct dates.
#' @return The input \code{data.table}, with a column \code{duration} added.
add_duration <- function(d) {
  td <- times_and_durations(d)
  # Can't test directly for POSIXct equality in current version of data.table.
  # Workaround from https://github.com/Rdatatable/data.table/issues/1008
  dur_from_time <- function(t) td[time >= t & time <= t, duration]
  d[, duration := dur_from_time(time), by = .(time)]
}

#' Creates a \code{data.table} with columns \code{time} and \code{duration}.
#' 
#' Creates a \code{data.table} with columns \code{time} and \code{duration},
#' in which the column \code{time} corresponds to the (unique) times in the 
#' input table. In the output, the most recent time given a duration of 1, 
#' the second most recent time given a duration of 2, and so on.
#' 
#' @param d A \code{data.table} containing at least a column \code{time},
#'        which is sortable. For example, could be POSIXct dates.
#' @return A new \code{data.table}, containing columns \code{time} (key column)
#'         and \code{duration}.
times_and_durations <- function(d) {
  times <- sort(unique(d[, time]), decreasing = FALSE)
  data.table(time = times, duration = rev(seq_along(times)), key = "time")
}

#' Calculate the log-likelihood ratios from given log-likelihoods.
#' 
#' From a \code{data.table} containing the log-likelihoods for all event types,
#' and for the null hypothesis of no event, calculates the log-likelihood ratios
#' by subtracting null log-likelihoods from the event log-likelihoods,
#' for each location, stream, and time.
#' 
#' @param loglikelihoods A \code{data.table} containing at least the columns 
#'        \code{event, location, stream, time} and \code{loglikelihood}.
#'        The first four columns must be key columns, in that order.
#'        The column \code{event} contains all event types (given e.g. as 
#'        integers or strings) and also the null hypothesis, as specified
#'        by the argument \code{null_name}.
#' @param null_name The identifier for the null hypothesis in the column
#'        \code{event} of the input argument \code{loglikelihoods}.
#'        E.g. \code{0L} if event types are specified as integers,
#'        or \code{"null"} if event types are specified as strings.
#' @return A \code{data.table} with key columns \code{location, event, stream,
#'         time}, and a column \code{llr} containing the log-likelihood ratios
#'         for each event type.
add_llr <- function(loglikelihoods, null_name) {
  input_keys <- c("event", "location", "stream", "time")
  if (any(getkeys(loglikelihoods)[1:4] != input_keys)) {
    stop("The key columns of the input table have to be ",
         "'event', 'location', 'stream', 'time'.")
  }
  keys <- c("location", "event", "stream", "time")
  loglikelihoods[
    event != null_name, 
    .(location = location, event = event, stream = stream, time = time, 
     llr = loglikelihood - loglikelihoods[event == null_name, loglikelihood])][,
    .SD, keyby = keys]
}

#' Get the set with the given index from an implicitly ordered set.
#' 
#' @param set_of_sets A \code{set} of \code{set}s. The elements of the outer
#'        set should be ordered, e.g. as it is when the elemets of the sets
#'        within it are integers.
#' @param index The index of the set you wish to be returned.
get_set <- function(set_of_sets, index) {
  i <- 1
  for (s in set_of_sets) {
    if (i == index) {
      return(s)
    }
    i <- i + 1
  }
}


#' Partition a set of regions.
#' 
#' Partition a set of regions such that each part contains about the same
#' number of locations, when the number of locations in each region for the
#' part are summed over all regions in the part.
#' @param regions A \code{set} of regions, each region itself being a \code{set}
#'        containing locations.
#' @param n_part An integer; the number of parts to split the \code{regions}
#'        into.
#' @return A list with two elements:
#'         \itemize{
#'         \item{partition} A list, each element of which is a \code{set} 
#'         containing one or more regions (\code{set} containing locations).
#'         \item{offsets} An integer vector containing offset numbers to the
#'         region numbering. For example, the first region in 
#'         \code{partition[i]} will have will be region number 
#'         \code{offset[i] + 1}.
#'         }
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
  total_n <- sum(n_locations)
  
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