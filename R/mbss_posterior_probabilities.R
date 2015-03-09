
# Warning: modifies argument by adding a column

#' Calculates the posterior log-probability of an event in a space window.
#' 
#' Calculates the posterior log-probability of an event of a given type \eqn{E_k}
#' for a given spatial region \eqn{S}\eqn{W},
#' for all event types and spatial regions.
#' Assumes a uniform distribution for the region prior probabilities.
#' This function \strong{modifies} the table passed as the first argument.
#' 
#' @param spacetime_llhs A \code{data.table}
#' @param event_logpriors A vector of the logarithms of the prior 
#'        probabilities of an event, ordered from event type 1 to K.
#' @param n_regions The number of regions.
#' @param data_logprob The logarithm of the marginal probability of the data.
spatial_logposterior <- function(spatial_llh,
                                 event_logpriors,
                                 n_regions,
                                 data_logprob) {
  spatial_llh[, posterior_logprob := llh + event_logpriors[event]
              - log(n_regions) - data_logprob]
  
}

#' Calculates the posterior log-probability of an event in a space-time window.
#' 
#' Calculates the posterior log-probability of an event of a given type \eqn{E_k}
#' for a given spatial region \eqn{S} and a given event duration \eqn{W},
#' for all event types, spatial regions and event durations up to \eqn{W_max}.
#' Assumes a uniform distribution for the region prior probabilities,
#' and likewise for the event durations.
#' This function \strong{modifies} the table passed as the first argument.
#' 
#' @param max_duration The number of time periods of the longest
#'        event duration considered.
#' @inheritParams spatial_llh
spacetime_logposterior <- function(spacetime_llh,
                                   event_logpriors,
                                   n_regions,
                                   max_duration,
                                   data_logprob) {
  spacetime_llh[, posterior_logprob := llh + event_logpriors[event]
                - log(n_regions) - log(max_duration) - data_logprob]
}

