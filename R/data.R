#' Longitude and latitude of New Mexico county seats.
#'
#' A dataset containing the longitude and latitude of the county seats of New
#' Mexico, except for Cibola county.
#'
#' @format A data table with 32 rows and 4 variables:
#' \describe{
#'   \item{county}{The name of the county. Of class character.}
#'   \item{seat}{The name of the county seat, i.e. the administrative center or
#'               seat of government. Of class factor with 32 levels.}
#'   \item{long}{The longitude of the county seat. Of class numeric.}
#'   \item{lat}{The latitude of the county seat. Of class numeric.}
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
#' @format A data table with 608 rows and 4 variables:
#' \describe{
#'   \item{year}{The year the cases were recorded. Of class integer.}
#'   \item{county}{The name of the county. Of class factor with 32 levels.}
#'   \item{population}{The population in that county and year. Of class 
#'                     numeric.}
#'   \item{count}{The number of brain cancer cases in that county and year. Of
#'                class integer.}
#' }
"NM_popcas" 

#' Data to plot the counties of New Mexico.
#' 
#' Map data for New Mexico. Was created using \code{ggplot2::map_data}.
#' @format A data table with 867 rows and 7 variables:
#' \describe{
#'   \item{long}{Longitude of county polygon corner. Of class numeric.}
#'   \item{lat}{Latitude of county polygon corner. Of class numeric.}
#'   \item{group}{Grouping by county. Of class numeric.}
#'   \item{order}{Order of the polygon corners. Of class integer.}
#'   \item{region}{Region is "new mexico" for all rows. Of class character.}
#'   \item{subregion}{The county name (with spaces). Of class character.}
#'   \item{county}{The county name (no spaces). Of class factor with 33 levels,
#'                 the county Cibola places last.}
#' }
"NM_map" 