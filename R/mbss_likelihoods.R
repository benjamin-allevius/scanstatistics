# Warning: these functions modify their arguments

#' Compute the log-likelihood log(\eqn{P(D|H_1(S,E_k))})
#'
#' Compute the log-likelihood log(\eqn{P(D|H_1(S,E_k))})
#' of the data given that an event of type \eqn{E_k} occurs in region 
#' \eqn{S}, for all such regions and for all event types,
#' considering events of time duration up to \eqn{W_{\max}}.
#' This function assumes that the event duration is uniformly distributed,
#' and that the event severity is uniformly distributed 
#' over its possible values.
#' 
#' @param lrsums A \code{data.table} with columns for region,
#'        location, event, and \code{llh},
#'        the log of the sum of likelihood ratios.
#' @param null_llh The log-likelihood for the complete data set under
#'        the null hypothesis of no event.
#' @param L The number of different possible event severities,
#'        or equivalently the number of possible impact factors.
#' @param max_duration The number of time periods of the longest
#'        event duration considered.
#' @return The input \code{data.table} with the \code{llh} column
#'         \strong{modified}, now containing the log-likelihoods
#'         log(\eqn{P(D|H_1(S,E_k))}).
#' @examples
#' lr_input <- data.table(region = rep(c(1,2,1,2), 2),
#'                        location = rep(c(1,2,3,3), 2), 
#'                        event = rep(1:2, each = 4),
#'                        llh = 1:8)
#' spatial_llh(lr_input, 
#'             null_llh = 1, 
#'             L = exp(2),
#'             max_duration = exp(-3))
spatial_llh <- function(full_llrs, null_llh, L, max_duration) {
  full_llrs[, 
            .(llh =  logsumexp(llr) + null_llh - log(L) - log(max_duration)),
            keyby = .(event, region)]
}

#' Compute log-likelihoods log(\eqn{P(D|H_1(S,E_k), W)}).
#'
#' Compute log-likelihoods log(P\eqn{(D|H_1(S,E_k), W)})
#' of the data given that an event of type \eqn{E_k} 
#' with time duration \eqn{W} occurs in spatial region \eqn{S},
#' for all durations \eqn{1 \le W \le W_{\max}},
#' all event types, and all regions.
#' Assumes the event severities are uniformly distributed over
#' \eqn{L} different values; same \eqn{L} for all event types.
#' 
#' @param full_llr A \code{data.table} with columns 
#'        \code{event, region, severity, time, llr}.
#'        The column \code{llr} contains the full log-likelihood ratios
#'        LR\eqn{_S^{k,l}(W)}.
#' @param null_llh The full log-likelihood (scalar value) 
#'        under the null hypothesis of no event. 
#' @param L The number of event severities, or equivalently the
#'        number of values the impact factor can take on.
#'        Assumed to be the same for all event types.
#' @return A \code{data.table} with columns \code{event, region, time, llh}.
#'         The column \code{llh} contains the log-likelihoods
#'         log(P\eqn{(D|H_1(S,E_k), W)}).
spacetime_llh <- function(full_llr, null_llh, L) {
  full_llr[, .(llh = logsumexp(llr) + null_llh - log(L)), 
           keyby = .(event, region, time)]
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
