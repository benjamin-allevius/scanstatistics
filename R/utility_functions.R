

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
