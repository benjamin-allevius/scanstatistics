

#' Computes the log-likelihood for each spatial region and event type.
#'
#' Computes the log-likelihood log(\eqn{P(D|H_1(S,E_k))})
#' of the data given that an event of type \eqn{E_k} occurs in region 
#' \eqn{S}, for all such regions and for all event types,
#' considering events of time duration up to \eqn{W_{\max}}.
#' This function assumes that the event duration is uniformly distributed,
#' and that the event severity is uniformly distributed 
#' over its possible values.
#' 
#' @inheritParams spacetime_llh_uniform
#' @param max_duration The number of time periods of the longest
#'        event duration considered.
#' @return A \code{data.table} with columns \code{region}, \code{event}, 
#'         \code{time}, and \code{llh}.
#'         The column \code{llh} contains the log-likelihoods
#'         for each region and event type.
spatial_llh_uniform <- function(full_llr, null_llh, L, max_duration) {
  full_llr[, .(llh =  logsumexp(llr) + null_llh - log(L) - log(max_duration)),
            keyby = .(region, event)]
}

#' Computes log-likelihood for each space-time region and event type.
#'
#' Computes the log-likelihoods log(P\eqn{(D|H_1(S,E_k), W)})
#' of the data given that an event of type \eqn{E_k} 
#' with time duration \eqn{W} occurs in spatial region \eqn{S},
#' for each duration \eqn{1 \le W \le W_{\max}},
#' each event type, and each region.
#' Assumes the event severities are uniformly distributed over
#' \eqn{L} different values; same \eqn{L} for all event types.
#' 
#' @param full_llr A \code{data.table} with columns 
#'        \code{region}, \code{event}, \code{severity}, \code{time}, 
#'        and \code{llr}. Keys should preferably be at 
#'        least \code{region, event}.        
#'        The column \code{llr} contains the full log-likelihood ratios.
#' @param null_llh The full log-likelihood (scalar value) 
#'        under the null hypothesis of no event. 
#' @param L The number of event severities, or equivalently the
#'        number of values the impact factor can take on.
#'        Assumed to be the same for all event types.
#' @return A \code{data.table} with columns \code{region}, \code{event}, 
#'         \code{time}, and \code{llh}.
#'         The column \code{llh} contains the log-likelihoods
#'         for each region, event type, and time.
spacetime_llh_uniform <- function(full_llr, null_llh, L) {
  full_llr[, .(llh = logsumexp(llr) + null_llh - log(L)), 
           keyby = .(region, event, time)]
}


#' Computes the full log-likelihood under the null hypothesis of no events.
#' 
#' Computes the log-likelihood log(P\eqn{(D|H_0}))
#' under the null hypotesis of no events taking place.
#' 
#' @param densities A \code{data.table} with column \code{density},
#'        containing the log-\emph{density} (log-pmf or log-pdf)
#'        for each observation in the data.
#'        Other columns should thus preferably be \code{location},
#'        \code{time}, and \code{stream}.
#' @return The log-likelihood log(\eqn{P(D|H_0})), a scalar value.
null_llh <- function(densities) {
  densities[, sum(density)]
}
