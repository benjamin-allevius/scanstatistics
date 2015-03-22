
# Note: we only calculate (log) likelihood ratios for time steps 0 <= t < W_max

# Input is a table with the log-likelihood ratio given the event effect, 
# and the log-probability of that effect,
# for each event, location, time, and stream.
# The function marginalizes over the effects.
# Calculates the log-likelihood ratios for each event, time, location,
# and stream.
event_llr <- function(llrs_given_effects) {
  llrs_given_effects[, .(llr = logsumexp(llr + effect_logprob)),
                     keyby = .(location, event, time, stream)]
}


#' Sums individual LLRs over all locations in each region.
#' 
#' For a given data stream \eqn{m}, time step \eqn{t}, region \eqn{S},
#' and event \eqn{E_k}, 
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
spacetime_llr <- function(event_llrs) {
  event_llrs[, .(llr = sum(llr)), by = .(region, event, time)][,
    .(time = time, llr = cumsum(llr)), by = .(region, event)]
} 


# eventtime_logprobs = matrix, rows representing event duration,
# columns representing event types,
# elements are log-probabilities of that duration for that event
spatial_llr <- function(st_llrs, eventtime_logprobs) {
  st_llrs[, .(llr = logsumexp(llr + eventtime_logprobs[time + 1, event])), 
          by = .(region, event)]
}


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