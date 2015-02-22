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
#' spatial_llh_uniform(lr_input, 
#'                     null_llh = 1, 
#'                     L = exp(2),
#'                     max_duration = exp(-3))
spatial_llh_uniform <- function(lrsums, 
                                null_llh, 
                                L,
                                max_duration) {
  lrsums[, llh := llh + null_llh - log(L) - log(max_duration)]
}

#' Compute log-likelihoods log(\eqn{P(D|H_1(S,E_k), W)})
#'
#' Compute log-likelihoods log(\eqn{P(D|H_1(S,E_k), W)})
#' of the data given that an event of type \eqn{E_k} 
#' with time duration \eqn{W} occurs in spatial region \eqn{S},
#' for all durations \eqn{1 \le W \le W_{\max}},
#' all event types, and all regions.
#' Assumes the event severities are uniformly distributed over
#' \eqn{L} different values; same \eqn{L} for all event types.
#' 
#' @param lrsums A \code{data.table} with columns for region,
#'        location, event, duration, and \code{llh},
#'        the log of the sum of likelihood ratios.
#' @param null_llh The log-likelihood of the complete data
#'        under the null hypothesis of no event. 
#' @param L The number of event severities, or equivalently the
#'        number of values the impact factor can take on.
#' @return The input \code{data.table} with the \code{llh} column
#'         \strong{modified}, now containing the log-likelihoods
#'         log(\eqn{P(D|H_1(S,E_k), W)}).
#' @examples
#' stlr_input <- data.table(region = rep(c(1,2,1,2), 4),
#'                          location = rep(c(1,2,3,3), 4),
#'                          event = rep(1:2, 2, each = 4),
#'                          duration = rep(1:2, each = 8),
#'                          llh = 1:16)
#' spacetime_llh_uniform(stlr_input, 
#'                       null_llh = 1, 
#'                       L = exp(2))
spacetime_llh_uniform <- function(lrsums, 
                                  null_llh, 
                                  L) {
  lrsums[, llh := llh + null_llh - log(L)]
}


#' Compute log-likelihood log(\eqn{P(D|H_0}))
#' 
#' Computes the log-likelihood log(\eqn{P(D|H_0}))
#' under the null hypotesis of no events.
#' 
#' @param llh_stream A \code{data.table} with \emph{key column} \code{stream}
#'        and column \code{llh_stream_value} which for each stream is the
#'        term(s) of the log-pmf or log-pdf that at most depends on the 
#'        particular stream \eqn{m}, not on the time index \eqn{t}
#'        nor the location \eqn{i}.
#' @param llh_rest A \code{data.table} with \emph{key column} \code{stream}
#'        and column \code{llh_stream_value} which for each stream is the
#'        sum over all time steps and locations, of the log-pmf or 
#'        log-pdf with the term(s) that depend at most on the stream \eqn{m}
#'        subtracted.
#' @return The log-likelihood log(\eqn{P(D|H_0})), a scalar value.
null_llh <- function(llh_stream, llh_rest) {
  llh_stream[llh_rest][, sum(llh_stream_value + llh_sum_st)]
}