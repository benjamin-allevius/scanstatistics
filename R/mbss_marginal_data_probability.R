




#' Log of the sum of spatial event likelihoods over all regions.
#' 
#' Computes the log of the sum of spatial event likelihoods 
#' P\eqn{(D|H_1(S,E_k))} over all regions \eqn{S}, for each event type \eqn{k}.
#' 
#' @param spatial_llhs A \code{data.table} of spatial event log-likelihoods,
#'        with columns \code{region}, \code{event}, and \code{llh}.
logsumexp_llh_over_regions <- function(spatial_llhs) {
  spatial_llhs[, .(llh = logsumexp(llh)), by = .(event)]
}