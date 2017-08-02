#' Longitude and latitude of New Mexico county seats.
#'
#' A dataset containing the longitude and latitude of the county seats of New
#' Mexico, except for Cibola county.
#'
#' @format A data frame with 33 rows and 7 variables:
#' \describe{
#'   \item{county}{Factor; the counties of New Mexico (no spaces).}
#'   \item{seat}{Character; the name of the county seat, i.e. the administrative 
#'               center or seat of government.}
#'   \item{area(km2)}{Numeric; the area in square kilometers of each county.}
#'   \item{seat_long}{Numeric; the longitude of the county seat.}
#'   \item{seat_lat}{Numeric; the latitude of the county seat.}
#'   \item{center_long}{Numeric; the longitude of the geographical center of
#'                      the county.}
#'   \item{center_lat}{Numeric; the latitude of the geographical center of the 
#'                     county.}
#' }
#' @source \url{https://en.wikipedia.org/wiki/List_of_counties_in_New_Mexico}
"NM_geo"

#' Population and brain cancer cases in New Mexico counties during 1973--1991.
#' 
#' A dataset containing the population count and number of brain cancer cases in 
#' the counties of New Mexico during the years 1973--1991. The population 
#' numbers are interpolations from the censuses conducted in 1973, 1982, and 
#' 1991. Interpolations were done using a quadratic function of time. Thus the
#' year-to-year changes are overly smooth but match the census numbers in the 
#' three years mentioned.
#' @format A data frame with 608 rows and 4 variables:
#' \describe{
#'   \item{year}{Integer; the year the cases were recorded.}
#'   \item{county}{Character; the name of the county (no spaces).}
#'   \item{population}{Integer; the population in that county and year.}
#'   \item{count}{Integer; the number of brain cancer cases in that county and 
#'                year.}
#' }
"NM_popcas" 

#' Data to plot the counties of New Mexico.
#' 
#' Map data for New Mexico. Was created using \code{ggplot2::map_data}.
#' @format A data frame with 867 rows and 7 variables:
#' \describe{
#'   \item{long}{Numeric; longitude of county polygon corner.}
#'   \item{lat}{Numeric; latitude of county polygon corner.}
#'   \item{group}{Numeric; grouping by county.}
#'   \item{order}{Numeric; order of the polygon corners.}
#'   \item{region}{Character; region is "new mexico" for all rows.}
#'   \item{subregion}{Character; the county name (with spaces).}
#'   \item{county}{Factor; the county name (no spaces).}
#' }
"NM_map" 
