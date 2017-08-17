
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
  
  if (is.null(population) && is.vector(counts)) {
    stop("Cannot reliably estimate baselines if counts are for a single time ",
         "point and no population data is available")
  }
  
  if (is.null(population) && is.matrix(counts) && nrow(counts) > 1) {
    return(tcrossprod(rowSums(counts), colSums(counts)) / total_count)
  }
  
  # If analysis is purely spatial:
  if (is.vector(counts) || nrow(counts) == 1) {
    if (is.matrix(population) && nrow(population) > 1) {
      stop("If counts is a vector, population should be too.")
    }
    return(matrix(population * total_count / sum(population), nrow = 1))
  }
  
  if (is.vector(counts)) {
    counts <- matrix(counts, nrow = 1)
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
#' Heuristically estimate the ZIP distribution Poisson mean parameters and the 
#' structural zero probabilities for each location and time point. Assumes the 
#' structural zero probability is constant over time for each location.
#' @param counts A matrix or vector of observed counts. Rows indicate time 
#'    (ordered from most recent) and columns indicate locations. If a vector,
#'    the elements are assumed to be the counts for each location. 
#' @param population A matrix or vector of populations for each location 
#'    (optional). If a matrix, should be of the same dimensions as 
#'    \code{counts}. If a vector, should be of the same length as the number of
#'    columns in \code{counts}.
#' @param min_p The minimum value you think possible for the structural zero
#'    probability.
#' @param min_mu The mimum value you think possible for the Poisson mean 
#'    parameter of the ZIP distribution (before adjusting for population size).
#' @return A list with two elements:
#'    \describe{
#'      \item{baselines}{A matrix of the same dimensions as \code{counts}.
#'                       If \code{counts} was a vector, a matrix with 1 row will
#'                       be returned.}
#'      \item{probs}{A matrix of the same dimensions as \code{counts}. If 
#'                   \code{counts} was a vector, a matrix with 1 row will be
#'                   returned.}
#'    }
#' @importFrom stats uniroot var
#' @keywords internal
estimate_zip_params <- function(counts, population = NULL, 
                                min_p = 0.001, min_mu = 0.3) {
  
  # Handle vector inputs
  if (is.vector(counts)) {
    counts <- matrix(counts, nrow = 1)
  }
  if (!is.null(population) && is.vector(population)) {
    population <- matrix(population, nrow = 1)
  }
  
  # If populations are given for a single timepoint, reshape to dim(counts)
  if (nrow(counts) > 1 && !is.null(population) && nrow(population) == 1) {
    population <- matrix(population, nrow(counts), ncol(counts), byrow = TRUE)
  }
  
  # ML estimation of ZIP parameters
  ml_est <- function(x) {
    # Function to solve for mu
    f <- function(mu) {
      prop_zeros <- max(min_p, sum(x == 0) / length(x))
      mean(x) * (1 - exp(-mu)) - mu * (1 - prop_zeros)
    } 
      
    par <- tryCatch( {
      # Try ML estimation
      mu_hat <- uniroot(f, interval = c(min_mu, max(min_mu * 2, max(x))))$root
      p_hat <- max(min_p, 1 - max(min_mu, mean(x)) / mu_hat)
      c(p_hat, mu_hat)
      
      # If ML fails, do method of moments
      }, error = function(cond) {
        a <- var(x) + mean(x)^2 - mean(x)
        b <- mean(x)
        d <- var(x) - mean(x)
        mu_hat <- ifelse(b == 0 || a < 0, min_mu, a / b)
        p_hat <- ifelse(d <= 0 || d / a <= 0, min_p, d / a)
        return(c(p_hat, mu_hat))
    })
    
    return(par)
  }
  
  
  if (nrow(counts) == 1) {
    par <- ml_est(counts)
    baselines <- matrix(par[2], nrow = 1, ncol = ncol(counts))
    probs <- matrix(par[1], nrow = 1, ncol = ncol(counts))
    # Adjust for population if applicable
    if (!is.null(population)) {
      baselines <- baselines * sum(counts) * population / sum(population)
    }
  } else {
    par <- apply(counts, 2, ml_est)
    baselines <- matrix(par[2, ], nrow = nrow(counts), ncol = ncol(counts),
                        byrow = TRUE)
    probs <- matrix(par[1, ], nrow = nrow(counts), ncol = ncol(counts),
                    byrow = TRUE)
    # Adjust for population if applicable
    if (!is.null(population)) {
      ff <- function(j) {
        baselines[, j] * sum(counts[, j]) * 
          population[, j] / sum(population[, j])
      }
      baselines <- sapply(seq_len(ncol(counts)), ff)
    }
  }
  
  list(baselines = baselines,probs = probs)
}
