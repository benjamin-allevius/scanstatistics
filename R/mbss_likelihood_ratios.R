
#' Computes log(\eqn{P(D | H_1(S,E_k), W, \theta_l^k) / P(D|H_0)})
full_lr <- function() {
  dt[time < maximum_duration,
     , 
     by = .(region, time, severity)]
}


#' Sums the location-dependent part of the LLR for all regions.
#' 
#' 
llr_region_sum <- function() {
  dt[time < maximum_duration, 
     , 
     by = .(region, event)]
}

#' Sums the location-independent part of the LLR for all regions.
#' 
#' 
llr_region_sum <- function() {
  dt[time < maximum_duration, 
     sum(likelihood) / number_of_severities, 
     by = .(region, event, time)]  
}