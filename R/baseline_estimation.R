
#' Estimate the \emph{baselines} (expected counts) by the Kulldorff method.
#' 
#' Estimates the baselines, which are the expected counts, by setting the 
#' expected count for a given time point and location to be the total count
#' for that time point multiplied by the proportion of all counts for that 
#' location.
#' 
#' @param count A \code{data.table} with columns 
#'        \code{stream, location, time, count}, keyed by the first three
#'        columns in that order.
#' @return A \code{data.table} with columns 
#'        \code{stream, location, time, count, baseline}. 
#'        Key columns are \code{stream, location, time} in that order.
kulldorff_baseline <- function(counts) {
  sum_by_time <- counts[, .(timesum = sum(count)), keyby = .(stream, location)]
  sum_by_loc <- counts[, .(locsum = sum(count)), keyby = .(stream, time)]
  sum_by_both <- counts[, .(totalsum = sum(count)), keyby = "stream"]
  
  timeprop <- merge(a, c, by = "stream")[, .(prop = timesum / totalsum), 
                                         keyby = .(stream, location)]
  baselines <- merge(b, p, by = c("stream"), allow.cartesian = T)[, 
    .(stream = stream, location = location, time = time, 
      baseline = locsum * prop)]
  setkeyv(baselines, c("stream", "location", "time"))
  merge(counts, baselines)
}