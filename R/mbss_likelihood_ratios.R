
# Note: we only calculate (log) likelihood ratios for time steps 0 <= t < W_max

#' Computes log(\eqn{P(D | H_1(S,E_k), W, \theta_l^k) / P(D|H_0)})
full_llr <- function() {
  dt[time < maximum_duration,
     , 
     by = .(region, time, severity)]
}


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


llr_time <- function() {
  DT5[, .(time = time, llr = cumsum(lq)), by = .(severity, event, region)]
}


sum_over_streams <- function(region_table, stream_table) {
  region_table[stream_table][, .
                             (lq = sum(lf, slg)),
                             keyby = .(time, severity, event, region)]
}


#' Sums the location-dependent part of the LLR for all regions.
#' 
#' 
sum_locations_in_region <- function(lg_table) {
  # sum_{i \in S} log g_{W-1}(i,m|x_{m,l}^k)
  lg_table[, .(slg = sum(lg)), 
           keyby = c("stream", "time", "severity", "event", "region")]
} 