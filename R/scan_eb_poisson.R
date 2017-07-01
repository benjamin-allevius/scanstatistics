

#' Calculate the expectation-based Poisson scan statistic.
#' @param counts A matrix of observed counts. Rows indicate time and are ordered
#'    from least recent (row 1) to most recent (row \code{nrow(counts)}).
#'    Columns indicate locations, numbered from 1 and up.
#' @param zones A list of integer vectors. Each vector corresponds to a single
#'    zone; its elements are the numbers of the locations in that zone.
#' @param baselines A matrix of the same dimensions as \code{counts}. Holds the
#'    Poisson mean parameter of the ZIP distribution for each observed count.
#'    Will be estimated if not supplied.
#' @param population A matrix or vector of populations for each location. Only
#'    needed if \code{baselines} are to be estimated and you want to account for
#'    the different populations in each location (and time).
#'    If a matrix, should be of the same dimensions as \code{counts}. If a
#'    vector, should be of the same length as the number of columns in
#'    \code{counts}.
#' @param n_mcsim A non-negative integer; the number of replicate scan
#'    statistics to generate in order to calculate a P-value.
#' @param max_only Boolean. If \code{FALSE} (default) the log-likelihood ratio
#'    statistic for each zone and duration is returned. If \code{TRUE}, only the
#'    largest such statistic (i.e. the scan statistic) is returned, along with
#'    the corresponding zone and duration.
#' @importFrom stats rpois
#' @importFrom ismev gum.fit
#' @importFrom reliaR pgumbel
#' @keywords internal
#' @examples
#' \dontrun{
#' set.seed(1)
#' # Create location coordinates, calculate nearest neighbors, and create zones
#' n_locs <- 50
#' max_duration <- 5
#' n_total <- n_locs * max_duration
#' geo <- matrix(rnorm(n_locs * 2), n_locs, 2)
#' knn_mat <- coords_to_knn(geo, 15)
#' zones <- knn_zones(knn_mat)
#'
#' # Simulate data
#' baselines <- matrix(rexp(n_total, 1/5), max_duration, n_locs)
#' counts <- matrix(rpois(n_total, as.vector(baselines)), max_duration, n_locs)
#'
#' # Inject outbreak/event/anomaly
#' ob_dur <- 3
#' ob_cols <- zones[[10]]
#' ob_rows <- max_duration + 1 - seq_len(ob_dur)
#' counts[ob_rows, ob_cols] <- matrix(
#'   rpois(ob_dur * length(ob_cols), 2 * baselines[ob_rows, ob_cols]), 
#'   length(ob_rows), length(ob_cols))
#' res <- scan_eb_poisson(counts = counts,
#'                        zones = zones,
#'                        baselines = baselines,
#'                        n_mcsim = 100,
#'                        max_only = FALSE,
#'                        rel_tol = 1e-3)
#' }
scan_eb_poisson <- function(counts,
                            zones,
                            baselines = NULL,
                            population = NULL,
                            n_mcsim = 0,
                            max_only = FALSE) {
  if (is.vector(counts)) {
    counts <- matrix(counts, nrow = 1)
  }
  if (!is.null(baselines) && is.vector(baselines)) {
    baselines <- matrix(baselines, nrow = 1)
  }
  

  # Estimate baselines and probs if not supplied
  if (is.null(baselines)) {
    baselines <- estimate_baselines(counts, population)
  } 
  
  # Reverse time order: most recent first
  counts <- counts[rev(seq_len(nrow(counts))), , drop = FALSE]
  baselines <- baselines[rev(seq_len(nrow(baselines))), , drop = FALSE]

  # Prepare zone arguments for C++
  zones_flat <- unlist(zones) - 1
  zone_lengths <- unlist(lapply(zones, length))
  num_locs <- ncol(counts)
  max_dur <- nrow(counts)
  num_zones <- length(zones)

  # Run analysis on observed counts
  scan <- scan_eb_poisson_cpp(counts, baselines,
                              zones_flat, zone_lengths,
                              num_locs, num_zones, max_dur, 
                              store_everything = !max_only)

  # Extract the most likely cluster (MLC)
  MLC <- scan[which.max(scan$score), ]

  # Make Monte Carlo replications of the scan statistic under the null hypothesis
  repl_stat <- numeric(n_mcsim)
  for (i in seq_len(n_mcsim)) {
    repl_stat[i] <- scan_eb_poisson_cpp(matrix(rpois(prod(dim(counts)),
                                                     as.vector(baselines)), 
                                               nrow(counts), ncol(counts)), 
                                        baselines,
                                        zones_flat, zone_lengths,
                                        num_locs, num_zones, max_dur, 
                                        store_everything = FALSE)$score
  }

  # Get P-values
  gumbel_pvalue <- NA
  MC_pvalue <- NA
  if (n_mcsim > 0) {
    gumbel_pvalue <- gumbel_pvalue(MLC$score, repl_stat, method = "ML")$pvalue
    MC_pvalue <- mc_pvalue(MLC$score, repl_stat)
  }

  list(MLC = list(zone_number = MLC$zone,
                  locations = zones[[MLC$zone]],
                  duration = MLC$duration,
                  score = MLC$score,
                  relative_risk = MLC$relrisk,
                  observed = counts[seq_len(MLC$duration),
                                    zones[[MLC$zone]]],
                  baselines = baselines[seq_len(MLC$duration),
                                        zones[[MLC$zone]]]),
       table = scan,
       replicate_statistics = repl_stat,
       MC_pvalue = MC_pvalue,
       Gumbel_pvalue = gumbel_pvalue)
}
