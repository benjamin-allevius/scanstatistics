
# Note: we only calculate (log) likelihood ratios for time steps 0 <= t < W_max




llr_sums_over_stream <- function() {
  DT[time < maximum_duration,
     , 
     by = .(region, time, severity)]
}


#' Sums the location-independent part of the LLR for all regions.
#' 
#' 
llr_region_sum <- function() {
  dt[time < maximum_duration, 
     sum(likelihood) / number_of_severities, 
     by = .(region, event, time)]  
}


#' Computes log(\eqn{P(D | H_1(S,E_k), W, \theta_l^k) / P(D|H_0)})
#' for all values 1 <= W <= W_max
full_llr <- function(lq_table) {
  lq_table[, .(time = time, llr = cumsum(lq)), 
           by = .(severity, event, region)]
}

#' Adds the region independent and dependent LLR parts together
#' over the data stream, then sums over the stream.
#' 
#' Both tables must have the following columns as the first 5 keys:
#' \code{c("stream", "time", "severity", "event")}.
#' 
#' @param region_table A \code{data.table} with
#' @param stream_table A \code{data.table} with
sum_over_streams <- function(region_table, 
                             stream_table) {
  region_table[stream_table][, .(lq = sum(lf, slg)),
                             keyby = .(time, severity, event, region)]
}


#' Sums the location-dependent part of the LLR for all regions.
#' 
#' @param lg_table A \code{data.table} containing columns
#'        \code{location}, \code{stream}, \code{time}, \code{severity}, 
#'        \code{event}, \code{region}, and \code{lg}.
#'        \code{lg} contains that part of the (full) log-likelihood
#'        ratio which depends on the locations contained in spatial region S.
sum_locations_in_region <- function(lg_table) {
  # sum_{i \in S} log g_{W-1}(i,m|x_{m,l}^k)
  lg_table[, .(slg = sum(lg)), 
           keyby = c("stream", "time", "severity", "event", "region")]
} 