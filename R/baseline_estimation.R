
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
#'        \code{stream, location, time, baseline}. Only \code{stream} will be
#'        a key column, but the table will implicity be ordered on
#'        \code{location} and \code{time} in the same way as the input table.
kulldorff_baseline <- function(counts) {
  merge(counts[, .(locsum = sum(count)), keyby = .(stream, time)], 
        merge(counts[, .(timesum = sum(count)), 
                     keyby = .(stream, location)], 
              counts[, .(totalsum = sum(count)), 
                     keyby = "stream"], 
              by = "stream")[, 
                             .(timeprop = timesum / totalsum), 
                             keyby = .(stream, location)], 
        by = c("stream"), allow.cartesian = T)[, 
          .(stream = stream, location = location, time = time, 
            baseline = locsum * timeprop)]
}