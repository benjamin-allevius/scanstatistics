
# Warning: some of these functions modify their arguments by adding a column

#' Calculates the posterior log-probability of 
#' each event type for each spatial region.
#' 
#' Calculates the posterior log-probability P(\eqn{H_1(S,E_k)|D}) of an 
#' event of a given type \eqn{E_k} for a given spatial region \eqn{S},
#' for all event types and spatial regions.
#' This function \strong{modifies} the table passed as the first argument,
#' by adding a column \code{posterior_logprob}.
#' @param spatial_llrs A \code{data.table} with columns \code{region, event}
#'    and \code{llr}, the latter containing the log-likelihood ratios for each 
#'    region and event type.
#' @param event_logpriors A vector of the logarithms of the prior probabilities 
#'    of an event, ordered from event type 1 to event type K (will be accessed 
#'    by index).
#' @param n_regions The total number of spatial regions.
#' @param data_logratio The logarithm of the ratio of the marginal probability 
#'    of data P(\eqn{D}) to the probability of the data under the null 
#'    hypothesis P(\eqn{D|H_0}).
spatial_logposterior <- function(spatial_llrs,
                                 event_logpriors,
                                 n_regions,
                                 data_logratio) {
  spatial_llrs[, posterior_logprob := llr + event_logpriors[event]
               - log(n_regions) - data_logratio]
  
}

#' Calculates the posterior log-probability of 
#' each event type for each spatial region and event duration.
#' 
#' Calculates the posterior log-probability of an 
#' event of a given type \eqn{E_k} for a given spatial 
#' region \eqn{S} and a given event duration \eqn{W},
#' for all event types, spatial regions and event durations considered.
#' This function \strong{modifies} the table passed as the first argument,
#' by adding a column \code{posterior_logprob}.
#' @inheritParams spatial_llr
#' @inheritParams spatial_logposterior
spacetime_logposterior <- function(spacetime_llrs,
                                   event_logpriors,
                                   dur_given_event_logprobs,
                                   n_regions,
                                   data_logratio) {
  spacetime_llrs[, posterior_logprob := llr + event_logpriors[event] + 
                   dur_given_event_logprobs[duration, event] - 
                   log(n_regions) - data_logratio,
                 by = .(region, event, duration)]
}


#' Log of the ratio [marginal probability of data] to 
#' [probability of data under null].
#' 
#' @param spatial_llrs A \code{data.table} with columns \code{region, event}
#'    and \code{llr}, the latter containing the log-likelihood ratios for each 
#'    region and event type.
#' @param event_logpriors A vector of the prior log-probabilities of each
#'    event type.
#' @param null_prior The prior probability of no event.
#' @param n_regions The total number of spatial regions.
#' @return The logarithm of the ratio of the marginal probability of data 
#'    P(\eqn{D}) to the probability of the data under the null hypothesis 
#'    P(\eqn{D|H_0}).
data_to_nulldata_logratio <- function(spatial_llrs, 
                                      event_logpriors,
                                      null_prior,
                                      n_regions) {
  log(null_prior + 
      exp(-log(n_regions) + 
           logsumexp(spatial_llrs[, .(rs = logsumexp(llr)), by = event][,
                                  event_logpriors[event] + rs])))
}

#' Computes the log of the marginal probability of the data.
#' 
#' Computes the logarithm of the marginal probability of the data.
#' @param data_logratio Output from \code{\link{data_to_nulldata_logratio}}.
#' @param null_llh The logarithm of the probability of the data under the
#'    null hypothesis of no events.
marginal_logprob_of_data <- function(data_logratio, null_llh) {
  data_logratio + null_llh
}

#' Calculates the posterior event probabilities.
#' 
#' @param dur_event_pjdist A \code{data.table} containing at least the columns
#'    \code{event} and \code{posterior_logprob}, the latter containing the log 
#'    of the posterior probability of the given event type (and given whatever 
#'    else is on the same row, in the other columns).
#' @return A \code{data.table} with columns \code{event} and 
#'    \code{event_posterior}, the latter containing the posterior probability 
#'    of the given event type.
posterior_event_probabilities <- function(dur_event_pjdist) {
  dur_event_pjdist[, .(event_posterior = sum(duration_event_posterior)), 
                   by = .(event)]
}

#' Calculates the posterior joint distribution of event type and duration.
#' 
#' @param spacetime_logposteriors A \code{data.table} containing at least the
#'    columns \code{event}, \code{duration}, and \code{posterior_logprob}.
#'    The latter contains the log of the posterior probability of the given 
#'    event type and duration, and given whatever else is on the same row, in 
#'    the remaining columns.
#' @return A \code{data.table} containing the columns \code{event},
#'    \code{duration}, and \code{duration_event_posterior}. The latter contains 
#'    the joint pmf of the event type and event duration.
posterior_duration_event_jdist <- function(spacetime_logposteriors) {
  spacetime_logposteriors[, 
    .(duration_event_posterior = exp(logsumexp(posterior_logprob))), 
    by = .(event, duration)]
}

#' Posterior conditional distribution of event duration given event type.
#' 
#' Given the posterior joint distribution of the event type and event duration,
#' calculates the posterior conditional distribution of the event duration
#' given the event type.
#' @param dur_event_pjdist A \code{data.table} containing the posterior
#'    joint distribution of the event duration and the event type.
posterior_duration_givn_event <- function(dur_event_pjdist) {
  dur_event_pjdist[, 
    .(duration = duration, 
      duration_condposterior = 
        duration_event_posterior / sum(duration_event_posterior)), 
    by = .(event)]
}