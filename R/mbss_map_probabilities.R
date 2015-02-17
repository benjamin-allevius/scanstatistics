
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
#' @param region_event_probabilities A \code{data.table},
#' containing the probabilities \eqn{P(H_1(S,E_k)|D)} for all 
#' regions \eqn{S\subseteq\{s_1,\ldots,s_N\}} and event 
#' types \eqn{E_k}, \eqn{k=1,\ldots,K}.
#' @return A \code{data.table} of probabilities \eqn{P(H_1(s_i,E_k)|D)} 
#' for each location \eqn{s_i} and event type \eqn{E_k}.
event_probability_map <- function(region_event_probabilities) {
    region_event_probabilities[, 
                               list(probability = sum(probability)), 
                               by = "event,location"]
}


#' Compute posterior probabilities \eqn{P(H_1(S,E_k)|D)}
#' and for each region \eqn{S\subseteq\{s_1,\ldots,s_N\}} and event
#' type \eqn{E_k}, \eqn{k=1,\ldots,K}.
#' 
#' @param region_event_loglikelihoods A \code{data.table} of 
#' @param region_event_logpriors A \code{data.table} of 
#' @param data_logprobability A scalar
#' @return The posterior probabilities \eqn{P(H_1(S,E_k)|D)}
#' for each region \eqn{S\subseteq\{s_1,\ldots,s_N\}} and event 
#' type \eqn{E_k}, \eqn{k=1,\ldots,K}.
#' @examples
#' reps <- data.table()
#' event_probability_map(reps)
region_event_probabilities <- function(region_event_loglikelihoods, 
                                       region_event_logpriors,
                                       data_logprobability) {
    
}