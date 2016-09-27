
#' Estimate the \emph{baselines} (expected counts) by the Kulldorff method.
#' 
#' Estimates the baselines, which are the expected counts, by setting the 
#' expected count for a given time point and location to be the total count
#' for that time point multiplied by the proportion of all counts for that 
#' location.
#' @param counts A \code{data.table} with columns 
#'    \code{stream, location, time, count}, keyed by the first three columns in 
#'    that order.
#' @return A \code{data.table} with columns 
#'    \code{stream, location, time, count, baseline}. Key columns are 
#'    \code{stream, location, time} in that order.
#' @keywords internal
kulldorff_baseline <- function(counts) {
  key_order <- c("stream", "location", "time")
  if (!all(getkeys(counts)[1:3] == key_order)) {
    stop("Key columns of input table must be 'stream', 'location', 'time'.")
  }
  sum_by_time <- counts[, .(timesum = sum(count)), keyby = .(stream, location)]
  sum_by_loc <- counts[, .(locsum = sum(count)), keyby = .(stream, time)]
  sum_by_both <- counts[, .(totalsum = sum(count)), keyby = "stream"]
  
  timeprop <- merge(sum_by_time, sum_by_both, by = "stream")[, 
                .(prop = timesum / totalsum), keyby = .(stream, location)]
  baselines <- merge(sum_by_loc, timeprop, 
                     by = c("stream"), allow.cartesian = TRUE)[, 
                 .(stream = stream, location = location, time = time, 
                   baseline = locsum * prop)]
  setkeyv(baselines, key_order)
  merge(counts, baselines)
}