# Note: we only calculate (log) likelihood ratios for durations 
# between 1 and W_max (inclusive)

#' Marginalizes log-likelihood ratio over event effects.
#' 
#' For a given location, duration, data stream, and event type affecing this 
#' \emph{unit of observation}, this function calculates the log-likelihood 
#' for the unit by marginalizing over the possible effects the event can have.
#' 
#' @param llrs_given_effects A \code{data.table} with columns 
#'        \code{location, event, duration, stream, llr, effect_logprob}.
#'        For a given row (i.e. location-duration-stream-event combination),
#'        the column \code{llr} contains the log-likelihood ratio
#'        corresponding to that unit of observation, assuming an event of
#'        that type occurs with a given effect. The column \code{effect_logprob}
#'        contains the log-probability of that specific event effect for that
#'        specific unit of observation.
#' @return A \code{data.table} with columns 
#'         \code{location, event, duration, stream, llr}; the first four being
#'         key columns (in that order) and the last containg the log-likelihood.
event_llr <- function(llrs_given_effects) {
  llrs_given_effects[, .(llr = logsumexp(llr + effect_logprob)),
                     keyby = .(location, event, duration, stream)]
}


#' Computes the log-likelihood ratio for all space-duration event combinations.
#' 
#' Computes the log-likelihood ratio for each spatial region \eqn{S},
#' event type \eqn{E_k}, and event duration \eqn{W}.
#' 
#' @param event_llrs A \code{data.table} with columns 
#'        \code{region, event, duration, stream, location}, and \code{llr}.
#'        The first 5 of these should preferably be key columns in the given
#'        order. The column \code{llr} contains the log-likelihood ratios.
#' @return A \code{data.table} with columns \code{region, event, duration, llr}.
spacetime_llr <- function(event_llrs) {
  event_llrs[, .(llr = sum(exp(llr))), by = .(region, event, duration)][,
    .(duration = duration, llr = log(cumsum(llr))), by = .(region, event)]
} 


#' Computes the log-likelihood ratio for all spatial regions and event types.
#' 
#' Computes the log-likelihood ratio for each spatial region \eqn{S} and
#' event type \eqn{E_k}.
#' 
#' @param spacetime_llrs A \code{data.table} with columns 
#'        \code{region, event, duration, llr}. The first three columns are 
#'        key columns in the given order, and column \code{llr} contains
#'        the space-duration log-likelihood ratios.
#' @param dur_given_event_logprobs A numeric matrix or data frame containing
#'        the conditonal log-probabilities of the duration of an event
#'        given the event type. Rows correspond to event durations and 
#'        columns correspond to events, so element \eqn{(i,j)} gives the
#'        log of the probability that the duration is \eqn{i},
#'        given that the event is of type \eqn{j}.
spatial_llr <- function(spacetime_llrs, dur_given_event_logprobs) {
  spacetime_llrs[, 
    .(llr = logsumexp(llr + dur_given_event_logprobs[duration, event])), 
    by = .(region, event)]
}
