
#' Compute posterior probability map \eqn{P(H_1(s_i)|D)}
#' for all spatial locations \eqn{s_i}.
#' 
#' @param event_probability_map A \code{data.table},
#' containing the probabilities \eqn{P(H_1(s_i,E_k)|D)} for each 
#' location \eqn{s_i} and event type \eqn{E_k}.
#' @return A vector of probabilities \eqn{P(H_1(s_i)|D)}.
#' @examples
#' epm <- data.table(event = rep(1:2, each = 5),
#' location = rep(1:5, 2),
#' probability = 1:10)
#' probability_map(epm)
probability_map <- function(event_probability_map) {
    event_probability_map[, 
                          list(probability = sum(probability)), 
                          by = location]
}


#' Compute probability map for all event types \eqn{P(H_1(s_i,E_k)|D)}
#' and all spatial locations \eqn{s_i}.
#' 
#' @param spatial_logposteriors A \code{data.table} containing the 
#' log-probabilities log(P\eqn{(H_1(S,E_k)|D)}) for all regions \eqn{S}
#' considered, and all event types \eqn{E_k}.
#' @param regions A \code{data.table} with key column \code{region} and
#'        second column \code{location}.
#' @return A \code{data.table} of probabilities \eqn{P(H_1(s_i,E_k)|D)} 
#' for each location \eqn{s_i} and event type \eqn{E_k}.
event_probability_map <- function(spatial_logposteriors, regions) {
  
}

