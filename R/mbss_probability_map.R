
#' Computes the posterior probability map for all locations.
#' 
#' Computes the posterior probability map \eqn{P(H_1(s_i)|D)}
#' for all locations \eqn{s_i}.
#' 
#' @param event_logprob_map A \code{data.table} with columns \code{location},
#'        \code{event}, and \code{posterior_logprob}.
probability_map <- function(event_logprob_map) {
  event_logprob_map[, .(probability = sum(exp(posterior_logprob))), 
                    by = .(location)]
}

#' Compute the probability map for all event types and spatial locations.
#' 
#' #' Compute the probability map for all event types \eqn{E_k}
#' and all locations \eqn{s_i}. The probability maps is,
#' for each such location and event type, the sum of the posterior probabilities
#' of that event type, over all regions \eqn{S} containing the location
#' \eqn{s_i}.
#' 
#' @inheritParams probability_map
event_probability_map <- function(event_logprob_map) {
  event_logprob_map[, .(probability = exp(posterior_logprob)), 
                    by = .(location, event)]
}

#' Compute the log-probability map for all event types and spatial locations.
#' 
#' Compute the logarithm of the probability map for all event types \eqn{E_k}
#' and all locations \eqn{s_i}. The probability maps is,
#' for each such location and event type, the sum of the posterior probabilities
#' of that event type, over all regions \eqn{S} containing the location
#' \eqn{s_i}.
#' 
#' @param spatial_logposteriors A \code{data.table} containing the key
#'        columns \code{region} and \code{event} (keyed in that order),
#'        and the column \code{posterior_logprob}, which contains the posterior
#'        log-probabilities log(P\eqn{(H_1(S,E_k)|D)}) for each region and 
#'        event. 
#' @param region_table A \code{data.table} with key column \code{region} and
#'        second column \code{location}.
#' @return A \code{data.table} of probabilities \eqn{P(H_1(s_i,E_k)|D)} 
#'         for each location \eqn{s_i} and event type \eqn{E_k}.
event_logprobability_map <- function(spatial_logposteriors, region_table) {
  region_table[spatial_logposteriors, allow.cartesian = TRUE][, 
    .(posterior_logprob = logsumexp(posterior_logprob)), 
    keyby = .(location, event)]
}

