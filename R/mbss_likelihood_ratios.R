

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