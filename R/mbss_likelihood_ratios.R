
# Note: we only calculate (log) likelihood ratios for time steps 0 <= t < W_max

#' Computes the full log-likelihood ratios.
#' 
#' Computes the full log-likelihood ratios; the log-likelihood ratio
#' for a given region, event, event severity, and time. 
#' 
#' @param llr_table A \code{data.table} containing the following columns 
#'        as keys in the given order: \code{region}, \code{event}, \code{time}, 
#'        and \code{severity}.
#'        Additionally, the column \code{llr} contains the sum of 
#'        log-likelihoods ratios over all data streams.
#' @return A \code{data.table} with key columns
#'        \code{region}, \code{event}, \code{severity}.
#'        The column \code{time} is not necessarily a key column,
#'        but is sorted from current time and backwards.
#'        The column \code{llr} contains the full log-likelihood ratios.
full_llr <- function(llr_table) {
  llr_table[, .(time = time, llr = cumsum(llr)), 
            by = .(region, event, severity)]
}

#' Sums the LLRs over all data streams.
#' 
#' For a given time step \eqn{t}, region \eqn{S},
#' and event \eqn{E_k} with severity \eqn{\theta_l^k}, 
#' sums the log-likelihood ratios over all data streams.
#' 
#' @param llr_table A \code{data.table} containing the following columns 
#'        as keys in the given order: \code{region}, \code{event}, \code{time}, 
#'        \code{severity}, \code{stream}.
#'        Additionally, the column \code{llr} contains the sum of 
#'        log-likelihood ratios over the locations in the given region.
#' @return A \code{data.table} with the same key columns as the input 
#'         \code{data.table}, less the column \code{stream}.
#'         The column \code{llr} now contains the sum
#'         of the column \code{llr} in the input \code{data.table}
#'         over all data streams.
sum_over_streams <- function(llr_table) {
  llr_table[, .(llr = sum(llr)), by = .(region, event, time, severity)]
}


#' Sums individual LLRs over all locations in each region.
#' 
#' For a given data stream \eqn{m}, time step \eqn{t}, region \eqn{S},
#' and event \eqn{E_k} with severity \eqn{\theta_l^k}, 
#' sums the log-likelihood ratios of the counts over the locations in \eqn{S}.
#' 
#' @param llr_table A \code{data.table} containing columns \code{location}
#'        and column \code{llr}, and the following columns as keys in
#'        the given order: \code{region}, \code{event}, \code{time}, 
#'        \code{severity}, \code{stream}.
#' @return A \code{data.table} with the same key columns as the input 
#'         \code{data.table}, with column \code{llr} containing the sum
#'         of the column \code{llr} in the input \code{data.table}
#'         over all locations in each region.
sum_locations_in_region <- function(llr_table) {
  llr_table[, .(llr = sum(llr)), by = .(region, event, time, severity, stream)]
} 