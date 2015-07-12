

#' Calculate the (logarithm of the) ZIP statistic for each space-time window.
#' 
#' Calculate the logarithm of the ZIP statistic for each space-time window.
#' @param table A \code{data.table} 
calculate_zip_statistic <- function(table, region) {
  table[, 
    statistic := sum(zip_stat_formula(relrisk, p, dstar, ddagger, mean, count)),
    by = .(region, duration)]
}


#' Formula for a term in the sum of the logarithm of the ZIP window statistic.
#' @param q The relative risk (scalar).
#' @param p Numeric vector of excess zero probabilities.
#' @param dstar Numeric vector of estimates of the excess zero indicators, under 
#'    the alternative hypothesis of an outbreak.
#' @param ddagger Numeric vector of estimates of the excess zero indicators, 
#'    under the null hypothesis of no outbreak.
#' @param mu Given/estimated Poisson means.
#' @param y Integer vector of observed counts.
zip_statistic_term <- function(q, p, dstar, ddagger, mu, y) {
  (dstar - ddagger) * (log(p) - log(1 - p) - y * log(mu) + lfactorial(y)) +
    (1 - dstar) * (y * log(q) - q * mu) + (1 - ddagger) * mu
}


#' Estimate the indicators of excess zeros for a ZIP distribution.
#' 
#' Given counts and (estimated) Poisson means and excess zero probabilities, 
#' this function estimates the indicator variable for an excess zero, for each
#' count.
#' @param p A numeric vector, each element being the given/estimated 
#'    probability of an excess zero for the corresponding count in \code{y}.
#' @param mu A numeric vector, each element being the given/estimated Poisson 
#'    meanfor the corresponding count in \code{y}. Of same length as \code{p}.
#' @param y An integer vector containing the observed counts. Of same length as 
#'    \code{p}.
#' @return A numeric vector, of same length as the input vectors.
estimate_d <- function(p, mu, y) {
  res <- rep(0, length(y))
  zero <- y == 0L
  res[zero] <- p[zero] / (p[zero] + (1 - p[zero]) * exp(-mu[zero]))
  res
}


#' Estimate the ZIP excess zero indicators under null hypothesis of no outbreak.
#' 
#' Estimate the excess zero indicator for a zero-inflated Poisson distribution,
#' under the null hypothesis of no outbreak. The estimate is added as a new 
#' column to the input table.
#' @param table A \code{data.table} with columns \code{p, mean, count}. \code{p}
#'    is the given/estimated probability of an excess zero, and \code{mean} is
#'    the estimated Poisson mean. \code{count} is the observed count.
#' @return The same table, \strong{modified} with a new column \code{ddagger}
estimate_d_dagger <- function(table) {
  table[, ddagger := estimate_d(p, mean, count)][]
}

#' Estimate the relative risk for an outbreak using a ZIP distribution.
#' 
#' Scalar estimate of the relative risk \eqn{q} for zero-inflated Poisson data.
#' @param d A vector of indicator variables for whether the corresponding count
#'    in the argument \code{y} is an excess zero or not. Can also be estimates 
#'    of these indicators.
#' @param mu A vector of given/estimated Poisson means, of same length as 
#'    \code{d}.
#' @param y An integer vector of the observed counts, of same length as 
#'    \code{d}.
#' @return A scalar, the estimated relative risk.
estimate_zip_relrisk <- function(d, mu, y) {
  max(1, sum(y * (1 - d)) / sum(mu * (1 - d)))
}

zip_relrisk_em <- function(p, mu, y, d_init = 0.5, tol = 0.01) {
  d_prev <- ifelse(y > 0, rep(0, length(y)), rep(d_init, length(y)))
  q <- estimate_zip_relrisk(d_prev, mu, y)
  d <- estimate_d(p, q * mu, y)
  
  left <- abs(d - d_prev) >= tol
  while (any(left)) {
    d_prev[left] <- d[left]
    q <- estimate_zip_relrisk(d, mu, y)
    d[left] <- estimate_d(p[left], q * mu[left], y[left])
    left[left] <- abs(d[left] - d_prev[left]) >= tol
  }
  # Return final q
  estimate_zip_relrisk(d, mu, y)
}