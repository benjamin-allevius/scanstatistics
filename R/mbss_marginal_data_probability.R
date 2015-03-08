

#' Log of event-related term of marginal probability of data.
#' 
#' Computes the log of the sum 
#' \deqn{N_S^{-1}\sum_{k=1}^K \sum_S P(D|H_1(S,E_k))},
#' which assumes a uniform prior distribution over regions.
#' 
#' @param event_llh A \code{data.table} with key column \code{event}
#'        and column \code{llh}, the latter containing the log of the
#'        sum of spatial event likelihoods over all regions.
#' @param event_logprior A \code{data.table} containing the logarithms
#'        of the prior event probabilities.
#' @param n_regions The total number of regions.
data_logprob_if_event <- function(event_llh, event_logprior, n_regions) {
  event_logprior[event_llh][, sum(event_logprob, llh)] - log(n_regions)
}


#' Log of the sum of spatial event likelihoods over all regions.
#' 
#' Computes the logarithm of the sum of spatial event likelihoods 
#' P\eqn{(D|H_1(S,E_k))} over all regions \eqn{S}, for each event type \eqn{k}.
#' 
#' @param spatial_llhs A \code{data.table} of spatial event log-likelihoods,
#'        with columns \code{region}, \code{event}, and \code{llh}.
logsumexp_llh_over_regions <- function(spatial_llhs) {
  spatial_llhs[, .(llh = logsumexp(llh)), keyby = .(event)]
}