# Functions in this file:
#   create_table

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
#' create_table(cols)
#' }
#' @keywords internal
#' @export
create_table <- function(col_list, keys = NULL) {
  data.table(do.call(expand.grid, c(col_list, stringsAsFactors = FALSE)),
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
#' @keywords internal
first_keys_match <- function(data_table, keys) {
  table_keys <- getkeys(data_table)
  if (is.null(table_keys) || length(keys) > length(table_keys)) {
    return(FALSE)
  } else {
    return(all(keys == table_keys[1:length(keys)]))
  }
}

#' Determine if table has all the given keys.
#' 
#' @inheritParams first_keys_match
#' @return \code{TRUE} if the table has all the given keys, \code{FALSE} 
#'    otherwise.
#' @keywords internal
haskeys <- function(data_table, keys) {
  all(keys %in% getkeys(data_table))
}

#' Get the keys from a data.table.
#' 
#' @param data_table A \code{data.table}.
#' @return NULL if the supplied \code{data.table} has no keys,
#'    else a character vector containing the keys.
#' @keywords internal
getkeys <- function(data_table) {
  attributes(data_table)$sorted
}

#' Extract the values of a \code{data.table} column by the column name.
#' 
#' This function extracts the values of a \code{data.table} column as a vector;
#' the column name is supplied as a character.
#' @param table A \code{data.table}.
#' @param colname The name of a column in the table.
#' @return The values of the column as a vector.
#' @keywords internal
get_column_values <- function(table, colname) {
  if (colname %notin% names(table)) {
    stop("The table does not contain a column named ", colname)
  }
  unname(table[, colname, with = FALSE][, unlist(.SD)])
}



#' Enumerate the unique, sorted values of a column, returned as a list.
#' 
#' Enumerate the unique and sorted values of a \code{data.table} column, and 
#' return the result as a list.
#' @param table A \code{data.table}.
#' @param colname The name of a column in the table. This column should be of
#'    class character or factor.
#' @return A list with the sorted and unique values of the column; the original
#'    values of the column are the list names, and the value for each name is 
#'    the number for that name (numbering starts at 1).
#' @keywords internal
enumerate_character <- function(table, colname) {
  if (colname %notin% names(table)) {
    stop(colname, " is not a name of a column in the table.")
  }
  if (sapply(table, class)[[colname]] %notin% c("character", "factor")) {
    warning("enumerate_character should only be used for character or factor ",
            "columns of a data.table.")
  }
  original <- get_column_values(table, colname) %>%
    unique %>%
    sort(decreasing = FALSE)
  enumerated <- seq_along(original)
  as.list(setNames(enumerated, original))
}

#' Get the name (as character) corresponding to the given number, from a list of
#' enumerated names.
#' 
#' @param enum_list A list as returned from \code{\link{enumerate_character}}.
#' @param number The number of the name you wish to extract.
#' @keywords internal
get_enumerated_character <- function(enum_list, number) {
  names(enum_list)[number]
}


#' Replaces a column with integers, its unique and sorted values enumerated.
#' 
#' Takes a \code{data.table} and replaces a single column (preferably) 
#' containing character or factor values with the enumeration of these values.
#' The enumeration is done according to the unique and sorted values of the 
#' original column elements. This function \strong{modifies} the input table.
#' @inheritParams enumerate_character
#' @return The input table \strong{modified}, with the given column replaced by 
#'    integer values as described.
#' @keywords internal
column_to_int <- function(table, colname) {
  timeclasses <- c("Date", "POSIXct", "POSIXlt")
  if (any(sapply(counts, class)[[colname]] %in% timeclasses)) {
    stop("Enumeration of time-related column classes not allowed.")
  }
  lookup <- enumerate_character(table, colname)
  get_int <- function(s) lookup[[s]]
  ints <- unname(vapply(table[, get(colname)], get_int, integer(1)))
  table[, c(colname) := list(ints)]
}
