

#' Create a \code{data.table} with all combinations of the supplied variables.
#' 
#' @param col_list A named list of the variables you want in your 
#'    \code{data.table}.
#' @param keys Character vector of one or more column names which is passed to 
#'    \code{\link[data.table]{setkey}}.
#' @return A \code{data.table} with all combinations of the variables supplied 
#'    in \code{col_list}.
#' @examples
#' \dontrun{
#' cols <- list(location = 1:2, time = 0:2, stream = 1:2)
#' table_creator(cols)
#' }
table_creator <- function(col_list, keys = NULL) {
  data.table(do.call(expand.grid, col_list),
             key = keys)
}

#' Do the first few keys of the supplied table match the supplied keys?
#' 
#' Checks if the first few key columns of the supplied \code{data.table} matches
#' those found in the supplied character vector, in the given order. Boolean 
#' output.
#' @param data_table A \code{data.table}.
#' @param keys A character vector.
#' @return \code{TRUE} if keys match, \code{FALSE} otherwise.
first_keys_match <- function(data_table, keys) {
  table_keys <- getkeys(data_table)
  if (is.null(table_keys) || length(keys) > length(table_keys)) {
    return(FALSE)
  } else {
    return(all(keys == table_keys[1:length(keys)]))
  }
}

haskeys <- function(data_table, keys) {
  all(keys %in% getkeys(data_table))
}

#' Get the keys from a data.table.
#' 
#' @param data_table A \code{data.table}.
#' @return NULL if the supplied \code{data.table} has no keys,
#'         else a character vector containing the keys.
getkeys <- function(data_table) {
  attributes(data_table)$sorted
}