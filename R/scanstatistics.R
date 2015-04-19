#' scanstatistics: A package for spatiotemporal event detection.
#'
#' The foo package provides three categories of important functions:
#' foo, bar and baz.
#' @section Foo functions:
#' The foo functions ...
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
  # Variables used unquoted inside functions
  "duration", 
  "W",
  "event",
  "location",
  "stream",
  "region",
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
  "regions_in_part",
  "totalsum",
  "loglikelihood",
  "variance",
  "aggregate_count",
  "aggregate_baseline",
  "score",
  "priority",
  "included_streams",
  "relative_risk",
  # data.table functions
  "data.table",
  "is.data.table",
  ".",
  ":=",
  "setkey",
  "setkeyv",
  ".SD",
  ".N",
  # foreach
  "%dopar%",
  "foreach"),
  package = "scanstatistics")