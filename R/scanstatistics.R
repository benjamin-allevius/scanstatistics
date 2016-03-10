#' scanstatistics: A package for spatiotemporal event detection.
#'
#' The scanstatistics package provides two categories of important functions:
#' data preparation functions, and the scan statistics themselves.
#' @section Data preparation functions:
#' These functions prepare your data for use. In particular, it helps you 
#' define the \emph{zones} which will be considered by the scan statistics.
#' @section Scan statistics:
#' These are the functions used for spacetime cluster detection.
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
  "duration", 
  "W",
  "event",
  "location",
  "stream",
  "zone",
  "posterior_prob",
  "posterior_logprob",
  "event_posterior",
  "llr",
  "rs",
  "count",
  "baseline",
  "overdispersion",
  "effect_logprob",
  "duration_condposterior",
  "num",
  "denom",
  "timesum",
  "locsum",
  "prop",
  "duration_event_posterior",
  "zones_in_part",
  "totalsum",
  "loglikelihood",
  "variance",
  "aggregate_count",
  "aggregate_baseline",
  # "score",
  "priority",
  "included_streams",
  "relative_risk",
  "phi",
  "relrisk",
  "statistic",
  "counts",
  "ddagger",
  "p",
  "score",
  "event_start_time",
## data.table functions---------------------------------------------------------
  "data.table",
  "is.data.table",
  ".",
  ":=",
  "setkey",
  "setkeyv",
  ".SD",
  ".N",
## foreach----------------------------------------------------------------------
  "%dopar%",
  "foreach"),
  package = "scanstatistics")