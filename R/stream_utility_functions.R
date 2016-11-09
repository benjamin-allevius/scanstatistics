# Functions in this file:
#   join_streamsets
#   create_streamset_table

#' Creates a new \code{data.table} from a table containing streams,
#' adding a column for sets of streams.
#' 
#' Takes a \code{data.table} containing the column \code{stream} and possibly 
#' other columns, and creates a new \code{data.table} with a column for 
#' \code{streamset} added to the columns in the supplied table, according to the 
#' argument \code{streamsets}. The key colums of the resulting \code{data.table} 
#' can be specified.
#' @param table A \code{data.table} with column \code{stream}
#'    and other columns (but none for \code{streamset}).
#' @param streamsets A list of stream sets, elements being vectors of integers
#'    representing different data streams.
#' @param keys Character vector of one or more column names; these columns
#'    are set as key columns in the output \code{data.table}.
#' @return A new \code{data.table} with a column for \code{streamset} added
#'    to the supplied table.
#' @examples
#' tab <- as.data.table(expand.grid(location = 1:2,
#'                                  duration = 1:3,
#'                                  stream = 1:2))
#' scanstatistics:::join_streamsets(tab, list(1, 2, 1:2), 
#'                                  keys = c("streamset", "duration"))
#' @keywords internal
join_streamsets <- function(table, streamsets, keys = c("streamsets")) {
  merge(x = table, 
        y = create_streamset_table(streamsets, keys = "stream"), 
        by = "stream", allow.cartesian = TRUE)[, .SD, keyby = keys]
}

#' Converts a list of stream sets to a \code{data.table} of stream sets and 
#' streams.
#' 
#' Supply a list of stream sets, with each element of the list being a vector
#' of streams. If the list is named, the output \code{data.table}
#' will have these names for the stream sets. Else, they will be labeled by 
#' integers from 1 to the length of the stream set list.
#' @param streamsets A list of stream sets, elements being vectors of integers
#'    representing different data streams.
#' @param keys Character vector of one or more column names which is passed 
#'    to \code{\link[data.table]{setkey}}.
#' @param offset An integer to offset the stream set numbering by, in case they
#'    are not named, and you want the numbering to start at \code{offset} + 1.
#' @examples 
#' \dontrun{
#' create_streamset_table(list(1L, 2L, 1:2))
#' create_streamset_table(sets::set(sets::set(1L), 
#'                      sets::set(2L), sets::as.set(1:2)))
#' create_streamset_table(list(1L, 2L, 1:2), keys = "stream")
#' create_streamset_table(list(1L, 2L, 1:2), keys = "streamset")
#' create_streamset_table(list(a = "x", b = "y", c = c("x", "y")))
#' }
#' @keywords internal
create_streamset_table <- function(streamsets, keys = NULL, offset = 0L) {
  streamset_names <- names(streamsets)
  if (is.null(streamset_names)) {
    streamset_names <- seq_along(streamsets) + offset
  }
  data.table(stream = unlist(streamsets, use.names = FALSE),
             streamset = rep(streamset_names, 
                             vapply(streamsets, length, integer(1))),
             key = keys)
}