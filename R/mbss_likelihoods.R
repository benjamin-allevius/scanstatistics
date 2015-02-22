
# Warning: these functions modify their arguments

#' Compute likelihoods \eqn{P(D|H_1(S,E_k))}
#' of the data given that an event of type \eqn{E_k} occurs in region 
#' \eqn{S\subseteq\{s_1,\ldots,s_N\}},
#' for all such regions and for all event types,
#' considering events of time duration up to \eqn{W_{\max}}.
#' 
#' @param region_event_duration_likelihoods A \code{data.table} of 
#' @param duration_probabilities A \code{data.table} of 
#' @return likelihoods \eqn{P(D|H_1(S,E_k))}
#' @examples
#' reps <- data.table()
#' event_probability_map(reps)
spatial_llh_uniform <- function(logsums, 
                                null_llh, 
                                n_severities,
                                max_duration) {
  logsums[, llh := llh + null_llh - log(n_severities) - log(max_duration)]
}


#' Compute log-likelihoods \eqn{\log P(D|H_1(S,E_k), W)}
#' of the data given that an event of type \eqn{E_k} 
#' with time duration \eqn{W} occurs in spatial region \eqn{S},
#' for all durations \eqn{1 \leq W \leq W_{\max}},
#' all event types, and all regions.
#' Assumes the event severities are uniformly distributed over
#' \eqn{L} different values; same \eqn{L} for all event types.
#' 
#' @param logsums A \code{data.table} of 
#' @param null_llh The log-likelihood of the complete data
#'        under the null hypothesis of no event. 
#' @param n_severities The number of event severities.
#' @return 
#' @examples
#' reps <- data.table()
#' event_probability_map(reps)
spacetime_llh_uniform <- function(logsums, 
                          null_llh, 
                          n_severities) {
  logsums[, llh := llh + null_llh - log(n_severities)]
}