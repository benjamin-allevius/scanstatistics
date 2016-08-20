
#' Extract the sorted names of all data streams from a \code{data.table}.
#' 
#' From a \code{data.table} with a column \code{stream}, extract the unique and
#' sorted values from this column as a vector. 
#' @param d A \code{data.table} containing a column \code{stream}.
#' @keywords internal
get_all_streams <- function(d) {
  if ("stream" %notin% names(d)) stop("data.table must have column 'stream'.")
  sort(unique(d[, stream]))
}

`%notin%` <- function(a, b) !(a %in% b)

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

#' Is the relative error between two numbers is less than the given tolerance?
#' 
#' Given two consecutive numbers in a sequence, return \code{TRUE} if the
#' relative change is positive but less than the given tolerance.
#' @param current A scalar; the most recent value of the sequence.
#' @param previous A scalar; the second most recent value of the sequence, or a
#'    reference value.
#' @param tol The tolerance, a positive scalar near zero.
#' @keywords internal
has_converged <- function(current, previous, tol = 0.01) {
  rel_change <- (current - previous) / abs(previous)
  rel_change > 0 && rel_change < tol
}



#' Add a column \code{duration} to a \code{data.table} with column \code{time}.
#' 
#' Adds a column \code{duration} to a \code{data.table} with column \code{time},
#' the most recent time given a duration of 1, the second most recent time 
#' given a duration of 2, and so on. This function \strong{modifies} the input
#' table.
#' @param d A \code{data.table} containing at least a column \code{time},
#'    which is sortable. For example, could be POSIXct dates.
#' @param keys A character vector for the columns you wish the output table to 
#'    be keyed by.
#' @param silent Should the table be modified and this function return nothing
#'    (silent = TRUE), or should the function return the modified table (silent
#'    = FALSE).
#' @return The input \code{data.table}, with a column \code{duration} added.
#' @keywords internal
add_duration <- function(d, keys = NULL, silent = TRUE) {
  td <- times_and_durations(d)
  # Can't test directly for POSIXct equality in current version of data.table.
  # Workaround from https://github.com/Rdatatable/data.table/issues/1008
  dur_from_time <- function(t) td[time >= t & time <= t, duration]
  d[, duration := dur_from_time(time), by = .(time)]
  if (!is.null(keys)) {
    setkeyv(d, keys)
  }
  if (!silent) {
    return(d)
  }
}

#' Creates a \code{data.table} with columns \code{time} and \code{duration}.
#' 
#' Creates a \code{data.table} with columns \code{time} and \code{duration},
#' in which the column \code{time} corresponds to the (unique) times in the 
#' input table. In the output, the most recent time given a duration of 1, the 
#' second most recent time given a duration of 2, and so on.
#' @param d A \code{data.table} containing at least a column \code{time}, which 
#'    is sortable. For example, could be POSIXct dates.
#' @return A new \code{data.table}, containing columns \code{time} (key column)
#'    and \code{duration}.
#' @keywords internal
times_and_durations <- function(d) {
  times <- sort(unique(d[, time]), decreasing = FALSE)
  data.table(time = times, duration = rev(seq_along(times)), key = "time")
}

#' Calculate the log-likelihood ratios from given log-likelihoods.
#' 
#' From a \code{data.table} containing the log-likelihoods for all event types,
#' and for the null hypothesis of no event, calculates the log-likelihood ratios
#' by subtracting null log-likelihoods from the event log-likelihoods,
#' for each location, stream, and event duration.
#' @param loglikelihoods A \code{data.table} containing at least the columns 
#'    \code{event, location, stream, duration} and \code{loglikelihood}. The 
#'    first four columns must be key columns, in that order. The column 
#'    \code{event} contains all event types (given e.g. as integers or strings) 
#'    and also the null hypothesis, as specified by the argument 
#'    \code{null_name}.
#' @param null_name The identifier for the null hypothesis in the column
#'    \code{event} of the input argument \code{loglikelihoods}. E.g. \code{0L} 
#'    if event types are specified as integers, or \code{"null"} if event types 
#'    are specified as strings.
#' @return A \code{data.table} with key columns \code{location, event, stream,
#'    duration}, and a column \code{llr} containing the log-likelihood ratios for 
#'    each event type.
#' @keywords internal
add_llr <- function(loglikelihoods, null_name) {
  input_keys <- c("event", "location", "stream", "duration")
  if (any(getkeys(loglikelihoods)[1:4] != input_keys)) {
    stop("The key columns of the input table have to be ",
         "'event', 'location', 'stream', 'duration'.")
  }
  keys <- c("location", "event", "stream", "duration")
  loglikelihoods[
    event != null_name, 
    .(location = location, event = event, stream = stream, duration = duration, 
     llr = loglikelihood - loglikelihoods[event == null_name, loglikelihood])][,
    .SD, keyby = keys]
}

#' Get the set with the given index from an implicitly ordered set.
#' 
#' @param set_of_sets A \code{set} of \code{set}s. The elements of the outer set 
#'    should be ordered, e.g. as it is when the elemets of the sets within it 
#'    are integers.
#' @param index The index of the set you wish to be returned.
#' @keywords internal
get_set <- function(set_of_sets, index) {
  i <- 1
  for (s in set_of_sets) {
    if (i == index) {
      return(s)
    }
    i <- i + 1
  }
}
