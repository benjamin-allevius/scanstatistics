
# Warning: modifies argument by adding a column

#' Calculates the posterior probability of an event in a space-time window.
#' 
#' Calculates the posterior probability of an event of a given type \eqn{E_k}
#' for a given spatial region \eqn{S} and a given event duration \eqn{W},
#' for all event types, spatial regions and event durations up to \eqn{W_max}.
#' This function \strong{modifies} the table passed as the first argument.
#' 
#' @param spacetime_llhs A \code{data.table}
#' @param event_logpriors A vector of the logarithms of the prior 
#'        probabilities of an event, ordered from event type 1 to K.
#' @param n_regions The number of regions.
#' @param max_duration The number of time periods of the longest
#'        event duration considered.
#' @param data_logprob The logarithm of the marginal probability of the data.
spacetime_logposterior <- function(spacetime_llh,
                                   event_logpriors,
                                   n_regions,
                                   max_duration,
                                   data_logprob) {
  spacetime_llh[, posterior_logprob := llh + event_logpriors[event]
                - log(n_regions) - log(max_duration) - data_logprob]
}

