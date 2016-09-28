#' scanstatistics: Space-time anomaly detection using scan statistics.
#'
#' The scanstatistics package provides two categories of important functions:
#' data preparation functions, and the scan statistics themselves.
#' @section Data preparation functions:
#' These functions prepare your data for use. In particular, it helps you 
#' define the \emph{zones} which will be considered by the scan statistics.
#' @section Scan statistics:
#' These are the functions used for space-time anomaly detection. Scan statistic
#' functions for univariate space-time data have a name that begins with 
#' \code{scan_} and functions for multivariate space-time data have a name that
#' begins with \code{mscan_}.
#' @docType package
#' @name scanstatistics
#' @import data.table
NULL



# To avoid problems with devtools::test()
# See https://github.com/hadley/devtools/issues/192
.datatable.aware=TRUE

# Hack based on Hadley Wickhams comment: 
# http://stackoverflow.com/a/12429344/897506
globalVariables(c(
## Variables used unquoted inside functions-------------------------------------
  "aggregate_baseline",
  "aggregate_count",
  "baseline",
  "count",
  "counts",
  "ddagger",
  "denom",
  "duration", 
  "duration_event_posterior",
  "duration_condposterior",
  "effect_logprob",
  "event",
  "event_posterior",
  "event_start_time",
  "included_streams",
  "llr",
  "location",
  "locsum",
  "loglikelihood",
  "mu",
  "num",
  "n_zones",
  "overdispersion",
  "p",
  "posterior_logprob",
  "posterior_prob",
  "priority",
  "prop",
  "relative_risk",
  "relative_score",
  "relrisk",
  "rs",
  "score",
  "statistic",
  "stream",
  "theta",
  "time",
  "timesum",
  "total_score",
  "totalsum",
  "variance",
  "W",
  "zone",  
  "zones_in_part",  
## data.table functions---------------------------------------------------------
  ".",
  ".N",
  ".SD",
  ":=",
  "data.table",
  "is.data.table",
  "setkey",
  "setkeyv",
# ## foreach----------------------------------------------------------------------
#   "%dopar%",
#   "foreach",
## misc ------------------------------------------------------------------------
  "offset",
  "setNames"),
  package = "scanstatistics",
  add = TRUE)