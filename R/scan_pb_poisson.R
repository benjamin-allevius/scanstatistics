#' Calculate the population-based Poisson scan statistic.
#' 
#' Calculate the population-based Poisson scan statistic devised by Kulldorff
#' (1997, 2001).
#' @param counts Either:
#'    \itemize{
#'      \item A matrix of observed counts. Rows indicate time and are ordered
#'            from least recent (row 1) to most recent (row 
#'            \code{nrow(counts)}). Columns indicate locations, numbered from 1 
#'            and up. If \code{counts} is a matrix, the optional argument 
#'            \code{population} should also be specified.
#'      \item A data frame with columns "time", "location", "count", 
#'            "population".
#'    }
#' @param zones A list of integer vectors. Each vector corresponds to a single
#'    zone; its elements are the numbers of the locations in that zone.
#' @param population Optional. A matrix or vector of populations for each 
#'    location and time point. Only needed if \code{baselines} are to be 
#'    estimated and you want to account for the different populations in each 
#'    location (and time). If a matrix, should be of the same dimensions as 
#'    \code{counts}. If a vector, should be of the same length as the number of 
#'    columns in \code{counts} (the number of locations).
#' @param n_mcsim A non-negative integer; the number of replicate scan
#'    statistics to generate in order to calculate a P-value.
#' @param max_only Boolean. If \code{FALSE} (default) the log-likelihood ratio
#'    statistic for each zone and duration is returned. If \code{TRUE}, only the
#'    largest such statistic (i.e. the scan statistic) is returned, along with
#'    the corresponding zone and duration.
#' @return A list which, in addition to the information about the type of scan
#'    statistic, has the following components:
#'    \describe{
#'      \item{MLC}{A list containing the number of the zone of the most likely
#'            cluster (MLC), the locations in that zone, the duration of the 
#'            MLC, the calculated score, and the relative risk inside and 
#'            outside the cluster. In order, the elements of this list are named  
#'            \code{zone_number, locations, duration, score, relrisk_in,
#'            relrisk_out}.}
#'      \item{observed}{A data frame containing, for each combination of zone 
#'            and duration investigated, the zone number, duration, score, 
#'            relative risks. The table is sorted by score with the top-scoring 
#'            location on top. If \code{max_only = TRUE}, only contains a single 
#'            row corresponding to the MLC.}
#'      \item{replicates}{A data frame of the Monte Carlo replicates of the scan 
#'            statistic (if any), and the corresponding zones and durations.}
#'      \item{MC_pvalue}{The Monte Carlo \eqn{P}-value.}
#'      \item{Gumbel_pvalue}{A \eqn{P}-value obtained by fitting a Gumbel 
#'            distribution to the replicate scan statistics.}
#'      \item{n_zones}{The number of zones scanned.}
#'      \item{n_locations}{The number of locations.}
#'      \item{max_duration}{The maximum duration considered.}
#'      \item{n_mcsim}{The number of Monte Carlo replicates made.}
#'    }
#' @references 
#'    Kulldorff, M. (1997). \emph{A spatial scan statistic}. Communications in 
#'    Statistics - Theory and Methods, 26, 1481–1496.
#'    
#'    Kulldorff, M. (2001). \emph{Prospective time periodic geographical disease 
#'    surveillance using a scan statistic}. Journal of the Royal Statistical 
#'    Society, Series A (Statistics in Society), 164, 61–72.
#' @importFrom dplyr arrange
#' @importFrom magrittr %<>%
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
                            population = NULL,
                            n_mcsim = 0,
                            max_only = FALSE) {
  if (is.data.frame(counts)) {
    # Validate input -----------------------------------------------------------
    if (any(c("time", "location", "count") %notin% names(counts))) {
      stop("Data frame counts must have columns time, location, count")
    }
    counts %<>% arrange(location, -time)
    # Create matrices ----------------------------------------------------------
    if ("population" %in% names(counts)) {
      population <- df_to_matrix(counts, "time", "location", "population")
    } 
    counts <- df_to_matrix(counts, "time", "location", "count")
  }
  
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
  args <- list(counts = counts, 
               baselines = baselines,
               zones = unlist(zones) - 1, 
               zone_lengths = unlist(lapply(zones, length)),
               store_everything = !max_only,
               num_mcsim = n_mcsim)
  
  # Run analysis on observed counts --------------------------------------------
  scan <- run_scan(scan_pb_poisson_cpp, args)
  
  MLC_row <- scan$observed[1, ]
  
  MLC_out <- list(zone_number = MLC_row$zone,
                  locations = zones[[MLC_row$zone]],
                  duration = MLC_row$duration,
                  score = MLC_row$score,
                  relrisk_in = MLC_row$relrisk_in,
                  relrisk_out = MLC_row$relrisk_out)
  
  structure(
    c(list(# General
      distribution = "Poisson",
      type = "population-based",
      setting = "univariate"),
      # MLC + analysis
      list(MLC = MLC_out),
      scan,
      # Configuration
      list(n_zones = length(zones),
           n_locations = ncol(counts),
           max_duration = nrow(counts),
           n_mcsim = n_mcsim)),
    class = "scanstatistic")
}
