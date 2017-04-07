

#' Calculate the population-based Poisson scan statistic.
#' @param counts A matrix of observed counts. Rows indicate time and are ordered
#'    from least recent (row 1) to most recent (row \code{nrow(counts)}).
#'    Columns indicate locations, numbered from 1 and up.
#' @param zones A list of integer vectors. Each vector corresponds to a single
#'    zone; its elements are the numbers of the locations in that zone.
#' @param baselines A matrix of the same dimensions as \code{counts}. Holds the
#'    Poisson mean parameter of the ZIP distribution for each observed count.
#'    Note that the sum of all elements in \code{baselines} must equal the sum
#'    of all elements in \code{counts}. If baselines is not supplied, it will be
#'    estimated.
#' @param population A matrix or vector of populations for each location. Only
#'    needed if \code{baselines} are to be estimated and you want to account for 
#'    the different populations in each location (and time).
#'    If a matrix, should be of the same dimensions as \code{counts}. If a
#'    vector, should be of the same length as the number of columns in
#'    \code{counts}.
#' @param n_mcsim A non-negative integer; the number of replicate scan
#'    statistics to generate in order to calculate a P-value.
#' @param gumbel Boolean. If \code{TRUE} and \code{n_mcsim > 0}, then a Gumbel
#'    distribution is fit to the replicate scan statistics and a \eqn{p}-value
#'    for the observed scan statistic is calculated using the fitted
#'    distribution. If \code{FALSE} and \code{n_mcsim > 0}, the empirical
#'    \eqn{p}-value is calculated instead, defined as 1 + the number of
#'    replicate scan statistics that were larger than the observed scan
#'    statistic, divided by \code{1 + n_mcsim}.
#' @param max_only Boolean. If \code{FALSE} (default) the log-likelihood ratio
#'    statistic for each zone and duration is returned. If \code{TRUE}, only the
#'    largest such statistic (i.e. the scan statistic) is returned, along with
#'    the corresponding zone and duration.
#' @importFrom stats rmultinom
#' @importFrom ismev gum.fit
#' @importFrom reliaR pgumbel
#' @keywords internal
#' @examples
#' \dontrun{
#' set.seed(1)
#' # Create location coordinates, calculate nearest neighbors, and create zones
#' geo <- matrix(rnorm(100), 50, 2)
#' knn_mat <- coords_to_knn(geo, 15)
#' zones <- knn_zones(knn_mat)
#'
#' # Simulate data
#' population <- 100 + matrix(2 * rt(50 * 5, df = 2), 5, 50)
#' baselines <- matrix(rexp(50 * 5, 1/5), 5, 50)
#' counts <- matrix(rpois(prod(dim(baselines)), as.vector(baselines)),
#'                  nrow(baselines), ncol(baselines))
#'
#' # Inject outbreak/event/anomaly
#' ob_dur <- 1:3
#' ob_zone <- zones[[10]]
#' counts[ob_dur, ob_zone] <- rpois(
#'   1, 2 * baselines[ob_dur, ob_zone])
#' res <- scan_pb_poisson(counts = counts,
#'                        population = population,
#'                        zones = zones,
#'                        n_mcsim = 100,
#'                        gumbel = TRUE,
#'                        max_only = FALSE)
#' }
scan_pb_poisson <- function(counts,
                            zones,
                            baselines = NULL,
                            population = NULL,
                            n_mcsim = 0,
                            gumbel = TRUE,
                            max_only = FALSE) {
  counts <- counts[rev(seq_len(nrow(counts))), ]
  total_count <- sum(counts)

  # Estimate baselines and probs if not supplied
  if (is.null(baselines)) {
    baselines <- estimate_baselines(counts, population)
  } else {
    baselines <- baselines[rev(seq_len(nrow(baselines))), ]
  }
  
  if (!isTRUE(all.equal(sum(baselines), total_count))) {
    stop("Sum of baselines must equal sum of counts.")
  }

  # Prepare zone arguments for C++
  zones_flat <- unlist(zones) - 1
  zone_lengths <- unlist(lapply(zones, length))

  # Run analysis on observed counts
  if (max_only) {
    scan <- scan_pb_poisson_cpp_max(counts, baselines,
                                    zones_flat, zone_lengths,
                                    total_count)
  } else {
    scan <- scan_pb_poisson_cpp(counts, baselines,
                                zones_flat, zone_lengths,
                                total_count)
  }

  # Extract the most likely cluster (MLC)
  MLC <- scan[which.max(scan$score), ]

  # Make Monte Carlo replications of the scan statistic under the null hypothesis
  repl_stat <- numeric(n_mcsim)
  for (i in seq_len(n_mcsim)) {
    repl_stat[i] <- scan_pb_poisson_cpp_max(
      matrix(rmultinom(1, total_count, as.vector(baselines)),
             nrow(baselines), ncol(baselines)),
      baselines,
      zones_flat, zone_lengths,
      total_count)$score
  }

  # Fit Gumbel distribution to Monte Carlo replicates
  gumbel_mu <- NA
  gumbel_sigma <- NA
  if (n_mcsim > 0 && gumbel) {
    gum_fit <- gum.fit(repl_stat, show = FALSE)
    gumbel_mu <- gum_fit$mle[1]
    gumbel_sigma <- gum_fit$mle[2]
  }

  # Get P-values if Monte Carlo simulations were made
  MC_pvalue <- NA
  gumbel_pvalue <- NA
  if (n_mcsim > 0) {
    MC_pvalue <- mc_pvalue(MLC$score, repl_stat)
    if (gumbel) {
      gumbel_pvalue <- pgumbel(MLC$score, gumbel_mu, gumbel_sigma,
                               lower.tail = FALSE)
    }
  }

  list(MLC = list(zone_number = MLC$zone,
                  locations = zones[[MLC$zone]],
                  duration = MLC$duration,
                  score = MLC$score,
                  relative_risk_in = MLC$relrisk_in,
                  relative_risk_out = MLC$relrisk_out,
                  observed = counts[seq_len(MLC$duration),
                                    zones[[MLC$zone]]],
                  baselines = baselines[seq_len(MLC$duration),
                                        zones[[MLC$zone]]]),
       table = scan,
       replicate_statistics = repl_stat,
       MC_pvalue = MC_pvalue,
       Gumbel_pvalue = gumbel_pvalue,
       Gumbel_mu = gumbel_mu,
       Gumbel_sigma = gumbel_sigma)
}
