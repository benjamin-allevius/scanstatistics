# Functions in this file:
#   mc_pvalue
#   print.scanstatistic
#   score_locations
#   top_clusters

#' Calculate the Monte Carlo \eqn{p}-value for a scan statistic.
#' 
#' Given an observed scan statistic \eqn{\lambda^*} and a vector of replicate 
#' scan statistics \eqn{\lambda_i}, \eqn{i=1,\ldots,R}, calculate the Monte 
#' Carlo \eqn{p}-value as
#' \deqn{
#'  \frac{1 + \sum_{i=1}^R \mathrm{I}(\lambda_i > \lambda^*)}{1 + R}
#' }
#' The function is vectorized, so multiple \eqn{p}-values can be calculated if
#' several scan statistics (e.g. statistics from secondary clusters) are 
#' supplied.
#' @param observed A scalar containing the observed value of the scan statistic,
#'    or a vector of observed values from secondary clusters.
#' @param replicates A vector of Monte Carlo replicates of the scan statistic.
#' @return The \eqn{p}-value or \eqn{p}-values corresponding to the observed 
#'    scan statistic(s).
#' @keywords internal
mc_pvalue <- function(observed, replicates) {
  if (length(replicates) == 0) {
    return(NULL)
  } else {
    f <- Vectorize(
      function(y) {
        (1 + sum(replicates > y)) / (1 + length(replicates))
        }
    )
    
    return(f(observed))
  }
}

#' Calculate the Gumbel \eqn{p}-value for a scan statistic.
#' 
#' Given an observed scan statistic \eqn{\lambda^*} and a vector of replicate 
#' scan statistics \eqn{\lambda_i}, \eqn{i=1,\ldots,R}, fit a Gumbel 
#' distribution to the replicates and calculate a \eqn{p}-value for the observed
#' statistic based on the fitted distribution.
#' \deqn{
#'  \frac{1 + \sum_{i=1}^R \mathrm{I}(\lambda_i > \lambda^*)}{1 + R}
#' }
#' The function is vectorized, so multiple \eqn{p}-values can be calculated if
#' several scan statistics (e.g. statistics from secondary clusters) are 
#' supplied.
#' @param observed A scalar containing the observed value of the scan statistic,
#'    or a vector of observed values from secondary clusters.
#' @param replicates A vector of Monte Carlo replicates of the scan statistic.
#' @param method Either "ML", for maximum likelihood, or "MoM", for method of 
#'    moments.
#' @return The \eqn{p}-value or \eqn{p}-values corresponding to the observed 
#'    scan statistic(s).
#' @importFrom ismev gum.fit
#' @importFrom reliaR pgumbel
#' @keywords internal
gumbel_pvalue <- function(observed, replicates, method = "ML") {
  # Fit Gumbel distribution to Monte Carlo replicates
  gumbel_mu <- NA
  gumbel_sigma <- NA
  if (method == "ML") {
    gum_fit <- gum.fit(replicates, show = FALSE)
    gumbel_mu <- gum_fit$mle[1]
    gumbel_sigma <- gum_fit$mle[2]
  } else {
    gumbel_sigma <- sqrt(6 * var(replicates) / pi^2)
    gumbel_mu <- mean(replicates) + digamma(1) * gumbel_sigma
  }
  
  pvalue <- pgumbel(observed, gumbel_mu, gumbel_sigma, lower.tail = FALSE)
  
  return(list(pvalue = pvalue, 
              gumbel_mu = gumbel_mu, 
              gumbel_sigma = gumbel_sigma))
}

#' Print a scanstatistic object.
#' 
#' Prints a scanstatistic object and returns it invisibly.
#' @param x A an object of class \code{scanstatistic}.
#' @param ... Further arguments passed to or from other methods.
#' @export
#' @keywords internal
print.scanstatistic <- function(x, ...) {
  if (x$type == "Bayesian") {
    cat(paste0(
      "Data distribution:                ", x$distribution, "\n",
      "Type of scan statistic:           ", x$type, "\n",
      "Setting:                          ", x$setting, "\n",
      "Number of locations considered:   ", x$n_locations, "\n",
      "Maximum duration considered:      ", x$max_duration, "\n",
      "Number of spatial zones:          ", x$n_zones, "\n",
      "Overall event probability:        ", x$posteriors$alt_posterior, "\n",
      "Probability of event in MLC:      ", round(x$MLC$posterior, 3), "\n",
      "Most likely event duration:       ", x$MLC$duration, "\n",
      "ID of locations in MLC:           ", toString(x$MLC$locations))
    )
  } else {
    cat(paste0(
      "Data distribution:                ", x$distribution, "\n",
      "Type of scan statistic:           ", x$type, "\n",
      "Setting:                          ", x$setting, "\n",
      "Number of locations considered:   ", x$n_locations, "\n",
      "Maximum duration considered:      ", x$max_duration, "\n",
      "Number of spatial zones:          ", x$n_zones, "\n",
      "Number of Monte Carlo replicates: ", x$n_mcsim, "\n",
      "Monte Carlo P-value:              ", ifelse(is.null(x$MC_pvalue), 
                                                   "NULL",
                                                   round(x$MC_pvalue, 3)), "\n",
      "Gumbel P-value:                   ", ifelse(is.null(x$Gumbel_pvalue), 
                                                   "NULL",
                                                 round(x$Gumbel_pvalue, 3)), "\n",
      "Most likely event duration:       ", x$MLC$duration, "\n",
      "ID of locations in MLC:           ", toString(x$MLC$locations))
      )
  }
  invisible(x)
}

#' Score each location over zones and duration.
#' 
#' For each location, compute the average of the statistic calculated for each
#' space-time window that the location is included in, i.e. average the 
#' statistic over both zones and the maximum duration.
#' @param x An object of class \code{scanstatistic}.
#' @param zones A list of integer vectors.
#' @return A \code{data.table} with the following columns:
#'    \describe{
#'      \item{location}{The locations (as integers).}
#'      \item{total_score}{For each location, the sum of all window statistics 
#'                         that the location appears in.}
#'      \item{n_zones}{The number of spatial zones that the location appears 
#'                     in.}
#'      \item{score}{The total score divided by the number of zones and the 
#'                   maximum duration.}
#'      \item{relative_score}{The score divided by the maximum score.}
#' }
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise
#' @importFrom tibble tibble
#' @export
#' @examples
#' \dontrun{
#' # Simple example
#' set.seed(1)
#' table <- data.frame(zone = 1:5, duration = 1, score = 5:1)
#' zones <- list(1:2, 1:3, 2:5, 4:5, c(1, 5))
#' x <- list(table = table, n_locations = 5, max_duration = 1, n_zones = 5)
#' score_locations(x, zones)
#' }
score_locations <- function(x, zones) {
  res <- tibble(location = seq_len(x$n_locations),
                score = 0,
                n_zones = 0)
  z_scores <- x$observed %>% 
    group_by(zone) %>% 
    summarise(score = sum(score)) %>%
    arrange(zone)
  
  if (nrow(z_scores) != length(zones)) stop("zones don't match x")
  
  for (i in seq_along(zones)) {
    for (location in zones[[i]]) {
      res[location, ]$score <- res[location, ]$score + z_scores$score[i]
      res[location, ]$n_zones <- res[location, ]$n_zones + 1
    }
  }
  res$score <- res$score / (x$n_zones * x$max_duration)
  res$relative_score <- res$score / max(res$score)
  return(res)
}

#' Get the top (non-overlappig) clusters.
#' 
#' Get the top \eqn{k} space-time clusters according to the statistic calculated
#' for each cluster (the maximum being the scan statistic). The default is to 
#' return the spatially non-overlapping clusters, i.e. those that do not have 
#' any locations in common.
#' @param x An object of class scanstatistics.
#' @param zones A list of integer vectors.
#' @param k An integer, the number of clusters to return.
#' @param overlapping Logical; should the top clusters be allowed to overlap in
#'    the spatial dimension? The default is \code{FALSE}.
#' @return A \code{tibble} with at most \eqn{k} rows, with columns 
#'    \code{zone, duration, score}. 
#' @export
#' @examples 
#' \dontrun{
#' set.seed(1)
#' table <- data.frame(zone = 1:5, duration = 1, score = 5:1)
#' zones <- list(1:2, 1:3, 2:5, c(1, 3), 4:5, c(1, 5))
#' top_clusters(list(table = table), zones, k = 4, overlapping = FALSE)
#' }
top_clusters <- function(x, zones, k = 5, overlapping = FALSE) {
  if (overlapping) {
    return(x$observed[seq_len(k), ])
  } else {
    row_idx <- c(1L, integer(k - 1))
    seen_locations <- zones[[1]]
    n_added <- 1L
    i <- 2L
    while (n_added < k && i <= nrow(x$observed)) {
      zone <- x$observed[i, ]$zone
      if (zone != x$observed[i-1, ]$zone && 
            length(intersect(seen_locations, zones[[zone]])) == 0) {
        seen_locations <- c(seen_locations, zones[[zone]])
        n_added <- n_added + 1L
        row_idx[n_added] <- i
      }
      i <- i + 1L
    }
    res <- x$observed[row_idx[row_idx > 0], ]
    res$MC_pvalue <- mc_pvalue(res$score, x$replicate_statistics$score)
    res$Gumbel_pvalue <- gumbel_pvalue(res$score, 
                                       x$replicates$score)$pvalue
    return(res)
  }
}
