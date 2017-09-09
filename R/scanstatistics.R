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
NULL

#' @useDynLib scanstatistics, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Hack based on Hadley Wickhams comment:
# http://stackoverflow.com/a/12429344/897506
globalVariables(c(
## Variables used unquoted inside functions-------------------------------------
  "location",
  "log_posterior",
  "score",
  "time",
  "zone"),
  package = "scanstatistics",
  add = TRUE)
