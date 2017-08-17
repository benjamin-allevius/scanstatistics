#' Calculate the expectation-based negative binomial scan statistic.
#'
#' Calculate the expectation-based negative binomial scan statistic devised by
#' Tango et al. (2011).
#' @param counts Either:
#'    \itemize{
#'      \item A matrix of observed counts. Rows indicate time and are ordered
#'            from least recent (row 1) to most recent (row 
#'            \code{nrow(counts)}). Columns indicate locations, numbered from 1 
#'            and up. If \code{counts} is a matrix, the optional matrix 
#'            arguments \code{baselines} and \code{thetas} should also be 
#'            specified.
#'      \item A data frame with columns "time", "location", "count", "baseline",
#'            "theta". See the description of the optional arguments 
#'            \code{baselines} and \code{thetas} below to see their definition.
#'    }
#' @param zones A list of integer vectors. Each vector corresponds to a single
#'    zone; its elements are the numbers of the locations in that zone.
#' @param baselines Optional. A matrix of the same dimensions as \code{counts}. 
#'    Holds the expected value parameter for each observed count. These 
#'    parameters are typically estimated from past data using e.g. GLM.
#' @param thetas Optional. A matrix of the same dimensions as \code{counts}, or 
#'    a scalar. Holds the dispersion parameter of the distribution, which is 
#'    such that if \eqn{\mu} is the expected value, the variance is 
#'    \eqn{\mu+\mu^2/\theta}. These parameters are typically estimated from past 
#'    data using e.g. GLM. If a scalar is supplied, the dispersion parameter is 
#'    assumed to be the same for all locations and time points.
#' @param type A string, either "hotspot" or "emerging". If "hotspot", the
#'    relative risk is assumed to be fixed over time. If "emerging", the
#'    relative risk is assumed to increase with the duration of the outbreak.
#' @param n_mcsim A non-negative integer; the number of replicate scan
#'    statistics to generate in order to calculate a \eqn{P}-value.
#' @param max_only Boolean. If \code{FALSE} (default) the statistic calculated
#'    for each zone and duration is returned. If \code{TRUE}, only the largest 
#'    such statistic (i.e. the scan statistic) is returned, along with the 
#'    corresponding zone and duration.
#' @return A list which, in addition to the information about the type of scan
#'    statistic, has the following components:
#'    \describe{
#'      \item{MLC}{A list containing the number of the zone of the most likely
#'            cluster (MLC), the locations in that zone, the duration of the 
#'            MLC, and the calculated score. In order, the 
#'            elements of this list are named  \code{zone_number, locations, 
#'            duration, score}.}
#'      \item{observed}{A data frame containing, for each combination of zone 
#'            and duration investigated, the zone number, duration, and score. 
#'            The table is sorted by score with the top-scoring location on top. 
#'            If \code{max_only = TRUE}, only contains a single row 
#'            corresponding to the MLC.}
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
#'    Tango, T., Takahashi, K. & Kohriyama, K. (2011), A space-time scan
#'    statistic for detecting emerging outbreaks, Biometrics 67(1), 106–115.
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
#'  baselines <- matrix(rexp(n_total, 1/5), max_duration, n_locs)
#'  thetas <- matrix(runif(n_total, 0.05, 3), max_duration, n_locs)
#'  counts <- matrix(rnbinom(n_total,  mu = baselines,  size = thetas), 
#'                   max_duration, n_locs)
#'
#' # Inject outbreak/event/anomaly
#' ob_dur <- 3
#' ob_cols <- zones[[10]]
#' ob_rows <- max_duration + 1 - seq_len(ob_dur)
#' counts[ob_rows, ob_cols] <- matrix(
#'   rnbinom(ob_dur * length(ob_cols), 
#'           mu = 2 * baselines[ob_rows, ob_cols],
#'           size = thetas[ob_rows, ob_cols]),
#'   length(ob_rows), length(ob_cols))
#' res <- scan_eb_negbin(counts = counts,
#'                       zones = zones,
#'                       baselines = baselines,
#'                       thetas = thetas,
#'                       type = "hotspot",
#'                       n_mcsim = 99,
#'                       max_only = FALSE)
#' }
scan_eb_negbin <- function(counts,
                           zones,
                           baselines,
                           thetas = 1,
                           type = c("hotspot", "emerging"),
                           n_mcsim = 0,
                           max_only = FALSE) {
  if (is.data.frame(counts)) {
    # Validate input -----------------------------------------------------------
    if (any(c("time", "location", "count", "baseline", "theta") %notin% 
              names(counts))) {
      stop("Data frame counts must have columns time, location, count, ",
           "baseline, theta.")
    }
    counts %<>% arrange(location, -time)
    # Create matrices ----------------------------------------------------------
    thetas <- df_to_matrix(counts, "time", "location", "theta")
    baselines <- df_to_matrix(counts, "time", "location", "baseline")
    counts <- df_to_matrix(counts, "time", "location", "count")
  }
  
  # Validate input -------------------------------------------------------------
  if (any(as.vector(counts) != as.integer(counts))) {
    stop("counts must be integer")
  }
  if (any(baselines <= 0)) stop("baselines must be positive")
  if (any(thetas <= 0)) stop("thetas must be positive")

  # Reshape arguments into matrices --------------------------------------------
  if (is.vector(counts)) {
    counts <- matrix(counts, nrow = 1)
  }
  if (!is.null(baselines) && is.vector(baselines)) {
    baselines <- matrix(baselines, nrow = 1)
  }
  
  if (is.vector(thetas)) {
    if (length(thetas) == 1) {
      thetas <- matrix(thetas, nrow(counts), ncol(counts))
    } else if (length(thetas) == ncol(counts)) {
      thetas <- matrix(thetas, nrow(counts), ncol(counts), byrow = TRUE)
    } else {
      stop("If thetas is supplied as a vector, it must be of the same length ",
           "as the number of locations.")
    }
  }

  # Reverse time order: most recent first --------------------------------------
  counts <- flipud(counts)
  baselines <- flipud(baselines)
  thetas <- flipud(thetas)

  # Prepare zone arguments for C++ ---------------------------------------------
  args <- list(counts = counts, 
               baselines = baselines,
               overdisp = 1 + baselines / thetas,
               zones = unlist(zones) - 1, 
               zone_lengths = unlist(lapply(zones, length)),
               store_everything = !max_only,
               num_mcsim = n_mcsim,
               score_hotspot = type[1] == "hotspot")

  # Run analysis on observed counts --------------------------------------------
  scan <- run_scan(scan_eb_negbin_cpp, args)
  
  MLC_row <- scan$observed[1, ]
  
  MLC_out <- list(zone_number = MLC_row$zone,
                  locations = zones[[MLC_row$zone]],
                  duration = MLC_row$duration,
                  score = MLC_row$score)

  structure(
    c(list(# General
      distribution = "negative binomial",
      type = "expectation-based",
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
