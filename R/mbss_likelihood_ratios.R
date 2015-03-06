
# Note: we only calculate (log) likelihood ratios for time steps 0 <= t < W_max

# NOT USED: can be used to improve speed if needed
#' Sum timesums of LR over severities
sum_lr_timesums_over_severities <- function(lr_timesums) {
  lr_timesums[, .(lr_sum = sum(lr_timesum)), keyby = .(event, region)]
}

# NOT USED: can be used to improve speed if needed
#' Sum LR over time
sum_lr_over_time <- function(lq_table) {
  lq_table[, .(lr_timesum = Reduce(function(x, y) x * (1 + y), 
                                    exp(lq), 
                                    right = TRUE)),
            keyby = .(event, region, severity)]
}


#' Computes the full log-likelihood ratios.
#' 
#' Computes the full log-likelihood ratios
#' LLR_S^{k,l}(W) = log(\eqn{P(D | H_1(S,E_k), W, \theta_l^k) / P(D|H_0)})
#' for all values 1 <= W <= W_max. 
#' 
#' @param lq_table A \code{data.table} 
full_llr <- function(lq_table) {
  lq_table[, .(time = time, llr = cumsum(lq)), 
           keyby = .(event, region, severity)]
}

#' Computes log-likelihood contribution for each time step.
#' 
#' Adds the region independent and dependent LLR parts together
#' over the data stream, then sums over the stream.
#' 
#' @section Important:
#' Both tables must have the following columns as the first 4 keys:
#' \code{c("stream", "time", "event", "severity")}.
#' 
#' @param region_table A \code{data.table} with
#' @param stream_table A \code{data.table} with
sum_over_streams <- function(region_table, 
                             stream_table) {
  region_table[stream_table][, .(lq = sum(lf, slg)),
                             keyby = .(time, event, severity, region)]
}


#' Sums the location-dependent part of the LLR for all regions.
#' 
#' @param lg_table A \code{data.table} containing columns
#'        \code{location}, \code{stream}, \code{time}, \code{severity}, 
#'        \code{event}, \code{region}, and \code{lg}.
#'        \code{lg} contains that part of the (full) log-likelihood
#'        ratio which depends on the locations contained in spatial region S.
#'        Table keys should preferably be
#'        c("stream", )
sum_locations_in_region <- function(lg_table) {
  lg_table[, .(slg = sum(lg)), 
           keyby = c("stream", "time", "event", "severity", "region")]
} 