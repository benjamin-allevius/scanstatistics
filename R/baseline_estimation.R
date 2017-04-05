
#' Estimate baselines based on observed counts.
#' 
#' Estimate the baselines (expected values) for the supplied counts.
#' @param counts A matrix of observed counts. Rows indicate time (ordered from 
#'    most recent) and columns indicate locations.
#' @param population A matrix or vector of populations for each location 
#'    (optional). If a matrix, should be of the same dimensions as 
#'    \code{counts}. If a vector, should be of the same length as the number of
#'    columns in \code{counts}.
#' @return A matrix of baselines of the same dimensions as \code{counts}.
#' @keywords internal
estimate_baselines <- function(counts, population = NULL) {
  total_count <- sum(counts)
  
  if (is.null(population)) {
    if (is.vector(counts)) {
      counts <- matrix(counts, nrow = 1)
    }
    return(tcrossprod(rowSums(counts), colSums(counts)) / total_count)
  }
  
  # If analysis is purely spatial:
  if (is.vector(counts)) {
    if (!is.vector(population)) {
      stop("If counts is a vector, population should be too.")
    }
    return(population * total_count / sum(population))
  }

  # If population is a vector, re-format to a matrix
  if (is.vector(population)) {
    population <- matrix(population,
                         nrow = nrow(counts),
                         ncol = length(population),
                         byrow = TRUE)
  }
  # Check that the same number of locations is implied by both arguments
  if (ncol(counts) != ncol(population)) {
    stop("The number of locations implied must be the same in the counts and ",
         "the population arguments.")
  }
  # Else return baselines as fractions of the total count
  timewise_pop <- rowSums(population)
  baselines <- matrix(0, nrow(counts), ncol(counts))
  for (i in seq_len(nrow(counts))) {
    baselines[i, ] <-
      population[i, ] * (total_count / (timewise_pop[i] * nrow(counts)))
  }
  return(baselines)
}

#' Estimate variances based on observed counts.
#' 
#' Estimate the variances for the supplied counts. It is assumed that variances
#' are constant over time for each location.
#' @param counts A matrix of observed counts. Rows indicate time (ordered from 
#'    most recent) and columns indicate locations.
#' @param baselines A matrix of the same dimensions as \code{counts} (optional).
#' @param population A matrix or vector of populations for each location 
#'    (optional). If a matrix, should be of the same dimensions as 
#'    \code{counts}. If a vector, should be of the same length as the number of
#'    columns in \code{counts}.
#' @param constant_dim An integer. If equal to 1, variances are assumed to be
#'    constant over time but different between locations. If equal to 2, 
#'    variances are assumed to vary over time but at each time point be equal 
#'    for all locations.
#' @return A matrix of variances of the same dimensions as \code{counts}.
#' @keywords internal
estimate_variances <- function(counts, 
                               baselines = NULL, 
                               population = NULL,
                               constant_dim = 1) {
  if (is.null(baselines)) {
    baselines <- estimate_baselines(counts, population)
  }
  
}

#' Estimate the parameters of a ZIP distribution.
#' 
#' Estimate the ZIP distribution Poisson mean parameters and the structural zero 
#' probabilities for each location and time point.
#' @param counts A matrix of observed counts. Rows indicate time (ordered from 
#'    most recent) and columns indicate locations.
#' @param population A matrix or vector of populations for each location 
#'    (optional). If a matrix, should be of the same dimensions as 
#'    \code{counts}. If a vector, should be of the same length as the number of
#'    columns in \code{counts}.
#' @param constant_dim An integer. If equal to 1, probabilities are assumed to 
#'    be constant over time but different between locations. If equal to 2, 
#'    probabilities are assumed to vary over time but at each time point be 
#'    equal for all locations.
#' @return A list with two elements:
#'    \describe{
#'      \item{baselines}{A matrix of the same dimensions as \code{counts}.}
#'      \item{probs}{A matrix of the same dimensions as \code{counts}.}
#'    }
#' @keywords internal
estimate_zip_params <- function(counts, population = NULL) {
  baselines <- estimate_baselines(counts, population)
  # Correct baselines somehow...
  probs <- baselines # fix
  
  list(baselines = baselines, probs = probs)
}


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