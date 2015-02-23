

# log-sum-exp trick
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
#' @example
#' cols <- list(location = 1:2, 
#'              time = 0:2, 
#'              stream = 1:2)
#' table_creator(cols)
table_creator <- function(col_list) {
  data.table(do.call(expand.grid, col_list))
}

# Takes a data.table of locations (locations_etc) and whatever other columns 
# (except regions), and creates a new table with a column for region added to
# the columns in locations_etc
# The key of locations_etc should be the column "location"
# @example
# cols <- table_creator(list(location = 1:2, time = 0:2, stream = 1:2))
# region_table_creator(cols, list(1, 2, 1:2))
region_table_creator <- function(locations_etc, regions) {
  
  region_names <- names(regions)
  if (is.null(region_names)) {
    region_names <- seq_along(regions)
  }
  
  DT <- data.table(location = unlist(regions, use.names = FALSE),
                   region = rep(region_names, 
                                vapply(regions, length, integer(1))))
  
  setkey(DT, "location")
  DT[locations_etc, allow.cartesian = TRUE]
}
