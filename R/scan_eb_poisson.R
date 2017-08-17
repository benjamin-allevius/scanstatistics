

#' Calculate the expectation-based Poisson scan statistic.
#' 
#' Calculate the expectation-based Poisson scan statistic devised by Neill et 
#' al. (2005).
#' @param counts Either:
#'    \itemize{
#'      \item A matrix of observed counts. Rows indicate time and are ordered
#'            from least recent (row 1) to most recent (row 
#'            \code{nrow(counts)}). Columns indicate locations, numbered from 1 
#'            and up. If \code{counts} is a matrix, the optional matrix argument 
#'            \code{baselines} should also be specified.
#'      \item A data frame with columns "time", "location", "count", "baseline". 
#'            Alternatively, the column "baseline" can be replaced by a column
#'            "population". The baselines are the expected values of the counts.
#'    }
#' @param zones A list of integer vectors. Each vector corresponds to a single
#'    zone; its elements are the numbers of the locations in that zone.
#' @param baselines Optional. A matrix of the same dimensions as \code{counts}. 
#'    Not needed if \code{counts} is a data frame. Holds the Poisson mean 
#'    parameter for each observed count. Will be estimated if not supplied 
#'    (requires the \code{population} argument). These parameters are typically 
#'    estimated from past data using e.g. Poisson (GLM) regression.
#' @param population Optional. A matrix or vector of populations for each 
#'    location. Not needed if \code{counts} is a data frame. If \code{counts} is
#'    a matrix, \code{population} is only needed if \code{baselines} are to be 
#'    estimated and you want to account for the different populations in each 
#'    location (and time). If a matrix, should be of the same dimensions as 
#'    \code{counts}. If a vector, should be of the same length as the number of 
#'    columns in \code{counts}.
#' @param n_mcsim A non-negative integer; the number of replicate scan
#'    statistics to generate in order to calculate a \eqn{P}-value.
#' @param max_only Boolean. If \code{FALSE} (default) the log-likelihood ratio
#'    statistic for each zone and duration is returned. If \code{TRUE}, only the
#'    largest such statistic (i.e. the scan statistic) is returned, along with
#'    the corresponding zone and duration.
#' @return A list which, in addition to the information about the type of scan
#'    statistic, has the following components:
#'    \describe{
#'      \item{MLC}{A list containing the number of the zone of the most likely
#'            cluster (MLC), the locations in that zone, the duration of the 
#'            MLC, the calculated score, and the relative risk. In order, the 
#'            elements of this list are named  \code{zone_number, locations, 
#'            duration, score, relative_risk}.}
#'      \item{observed}{A data frame containing, for each combination of zone 
#'            and duration investigated, the zone number, duration, score, 
#'            relative risk. The table is sorted by score with the top-scoring 
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
#'    Neill, D. B., Moore, A. W., Sabhnani, M. and Daniel, K. (2005). 
#'    \emph{Detection of emerging space-time clusters}. Proceeding of the 
#'    eleventh ACM SIGKDD international conference on Knowledge discovery in 
#'    data mining - KDD ’05, 218.
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
#'                        n_mcsim = 99,
#'                        max_only = FALSE)
#' }
scan_eb_poisson <- function(counts,
                            zones,
                            baselines = NULL,
                            population = NULL,
                            n_mcsim = 0,
                            max_only = FALSE) {
  if (is.data.frame(counts)) {
    # Validate input -----------------------------------------------------------
    if (any(c("time", "location", "count") %notin% names(counts))) {
      stop("Data frame counts must have columns time, location, count")
    }
    if (all(c("baseline", "population") %notin% names(counts))) {
      stop("For data frame counts, either population or baseline must be ",
           "specified.")
    }
    counts %<>% arrange(location, -time)
    # Create matrices ----------------------------------------------------------
    if ("baseline" %notin% names(counts)) {
      baselines <- NULL
      population <- df_to_matrix(counts, "time", "location", "population")
    } else {
      baselines <- df_to_matrix(counts, "time", "location", "baseline")
    }
    counts <- df_to_matrix(counts, "time", "location", "count")
  }
  
  # Validate input -------------------------------------------------------------
  if (any(as.vector(counts) != as.integer(counts))) {
    stop("counts must be integer")
  }
  if (!is.null(baselines) && any(baselines <= 0)) {
    stop("baselines must be positive")
  }
  
  # Reshape into matrices ------------------------------------------------------
  if (is.vector(counts)) {
    counts <- matrix(counts, nrow = 1)
  }
  if (!is.null(baselines) && is.vector(baselines)) {
    baselines <- matrix(baselines, nrow = 1)
  }

  # Estimate baselines if not supplied -----------------------------------------
  if (is.null(baselines)) {
    baselines <- estimate_baselines(counts, population)
  } 
  
  # Reverse time order: most recent first --------------------------------------
  counts <- flipud(counts)
  baselines <- flipud(baselines)

  # Prepare zone arguments for C++ ---------------------------------------------
  args <- list(counts = counts, 
               baselines = baselines,
               zones = unlist(zones) - 1, 
               zone_lengths = unlist(lapply(zones, length)),
               store_everything = !max_only,
               num_mcsim = n_mcsim)


  # Run analysis on observed counts --------------------------------------------
  scan <- run_scan(scan_eb_poisson_cpp, args)
  
  MLC_row <- scan$observed[1, ]
  
  MLC_out <- list(zone_number = MLC_row$zone,
                  locations = zones[[MLC_row$zone]],
                  duration = MLC_row$duration,
                  score = MLC_row$score,
                  relative_risk = MLC_row$relrisk)

  structure(
    c(list(# General
           distribution = "Poisson",
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
