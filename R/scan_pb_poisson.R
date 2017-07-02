#' Calculate the population-based Poisson scan statistic.
#' 
#' Calculate the population-based Poisson scan statistic devised by Kulldorff
#' (1997, 2001).
#' @param counts A matrix of observed counts. Rows indicate time and are ordered
#'    from least recent (row 1) to most recent (row \code{nrow(counts)}).
#'    Columns indicate locations, numbered from 1 and up.
#' @param zones A list of integer vectors. Each vector corresponds to a single
#'    zone; its elements are the numbers of the locations in that zone.
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
#' @return A list with the following components:
#'    \describe{
#'      \item{MLC}{A list containing the number of the zone of the most likely
#'            cluster (MLC), the locations in that zone, the duration of the 
#'            MLC, the calculated score, the relative risk inside and outside
#'            the MLC, and matrices of the observed counts, population, and
#'            estimated baselines for each location and time point in the MLC.}
#'      \item{table}{A data frame containing, for each combination of zone and
#'            duration investigated, the zone number, duration, score, relative 
#'            risk. If \code{max_only = TRUE}, only contains a single row 
#'            corresponding to the MLC.}
#'      \item{replicate_statistics}{A vector of the Monte Carlo replicates of
#'            the scan statistic, if any (otherwise empty).}
#'      \item{MC_pvalue}{The Monte Carlo \eqn{P}-value.}
#'      \item{Gumbel_pvalue}{A \eqn{P}-value obtained by fitting a Gumbel 
#'            distribution to the replicate scan statistics.}
#'      \item{n_zones}{The number of zones scanned.}
#'      \item{n_locations}{The number of locations.}
#'      \item{max_duration}{The maximum duration considered.}
#'    }
#' @references 
#'    Kulldorff, M. (1997). \emph{A spatial scan statistic}. Communications in 
#'    Statistics - Theory and Methods, 26, 1481–1496.
#'    
#'    Kulldorff, M. (2001). \emph{Prospective time periodic geographical disease 
#'    surveillance using a scan statistic}. Journal of the Royal Statistical 
#'    Society, Series A (Statistics in Society), 164, 61–72.
#' @importFrom stats rmultinom
#' @importFrom ismev gum.fit
#' @importFrom reliaR pgumbel
#' @export
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
#' population <- matrix(rnorm(n_total, 100, 10), max_duration, n_locs)
#' counts <- matrix(rpois(n_total, as.vector(population) / 20), 
#'                  max_duration, n_locs)
#'
#' # Inject outbreak/event/anomaly
#' ob_dur <- 3
#' ob_cols <- zones[[10]]
#' ob_rows <- max_duration + 1 - seq_len(ob_dur)
#' counts[ob_rows, ob_cols] <- matrix(
#'   rpois(ob_dur * length(ob_cols), 2 * population[ob_rows, ob_cols] / 20), 
#'   length(ob_rows), length(ob_cols))
#' res <- scan_pb_poisson(counts = counts,
#'                        zones = zones,
#'                        population = population,
#'                        n_mcsim = 99,
#'                        max_only = FALSE)
#' }
scan_pb_poisson <- function(counts,
                            zones,
                            population,
                            n_mcsim = 0,
                            max_only = FALSE) {
  # Validate input -------------------------------------------------------------
  if (any(as.vector(counts) != as.integer(counts))) {
    stop("counts must be integer")
  }
  if (any(population <= 0)) stop("population must be positive")
  
  # Reshape into matrices ------------------------------------------------------
  if (is.vector(counts)) {
    counts <- matrix(counts, nrow = 1)
  }
  if (is.vector(population)) {
    population <- matrix(population, nrow = 1)
  }
  
  # Estimate baselines ---------------------------------------------------------
  baselines <- estimate_baselines(counts, population)
  
  # Reverse time order: most recent first --------------------------------------
  counts <- flipud(counts)
  population <- flipud(population)
  baselines <- flipud(baselines)
  
  # Prepare zone arguments for C++ ---------------------------------------------
  zones_flat <- unlist(zones) - 1
  zone_lengths <- unlist(lapply(zones, length))
  num_locs <- ncol(counts)
  max_dur <- nrow(counts)
  num_zones <- length(zones)
  total_count <- sum(counts)
  
  # Run analysis on observed counts --------------------------------------------
  scan <- scan_pb_poisson_cpp(counts, baselines, total_count,
                              zones_flat, zone_lengths,
                              num_locs, num_zones, max_dur, 
                              store_everything = !max_only)
  
  # Extract the most likely cluster (MLC)
  MLC <- scan[which.max(scan$score), ]
  
  # Make MC replications of the scan statistic under the null hypothesis
  repl_stat <- numeric(n_mcsim)
  for (i in seq_len(n_mcsim)) {
    repl_stat[i] <- scan_pb_poisson_cpp(
      matrix(rmultinom(1, total_count, as.vector(baselines)), 
             nrow(baselines), ncol(baselines)), 
      baselines, total_count,
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
  
  MLC_counts <- counts[seq_len(MLC$duration), zones[[MLC$zone]], drop = FALSE]
  MLC_basel <- baselines[seq_len(MLC$duration), zones[[MLC$zone]], drop = FALSE]
  MLC_pop <- population[seq_len(MLC$duration), zones[[MLC$zone]], drop = FALSE]
  
  list(MLC = list(zone_number = MLC$zone,
                  locations = zones[[MLC$zone]],
                  duration = MLC$duration,
                  score = MLC$score,
                  relative_risk_in = MLC$relrisk_in,
                  relative_risk_out = MLC$relrisk_out,
                  observed = flipud(MLC_counts),
                  population = flipud(MLC_pop),
                  baselines = flipud(MLC_basel)),
       table = scan,
       replicate_statistics = repl_stat,
       MC_pvalue = MC_pvalue,
       Gumbel_pvalue = gumbel_pvalue,
       n_zones = length(zones),
       n_locations = ncol(counts),
       max_duration = nrow(counts))
}
