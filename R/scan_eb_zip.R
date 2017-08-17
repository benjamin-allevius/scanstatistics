

#' Calculate the expectation-based ZIP scan statistic.
#' 
#' Calculates the expectation-based scan statistic. See details below.
#' @param counts Either:
#'    \itemize{
#'      \item A matrix of observed counts. Rows indicate time and are ordered
#'            from least recent (row 1) to most recent (row 
#'            \code{nrow(counts)}). Columns indicate locations, numbered from 1 
#'            and up. If \code{counts} is a matrix, the optional matrix 
#'            arguments \code{baselines} and \code{probs} should also be 
#'            specified.
#'      \item A data frame with columns "time", "location", "count", "baseline",
#'            "prob". The baselines are the expected values of the counts, and 
#'            "prob" are the structural zero probabilities of the counts. If
#'            "baseline" and "prob" are not found as columns, their values are
#'            estimated in a \emph{very} heuristic fashion (not recommended).
#'            If population numbers are available, they can be included in a 
#'            column "population" to help with the estimation.
#'    }
#' @param zones A list of integer vectors. Each vector corresponds to a single
#'    zone; its elements are the numbers of the locations in that zone.
#' @param baselines Optional. A matrix of the same dimensions as \code{counts}. 
#'    Holds the Poisson mean parameter of the ZIP distribution for each observed 
#'    count. These parameters are typically estimated from past data using e.g. 
#'    ZIP regression.
#' @param probs Optional. A matrix of the same dimensions as \code{counts}. 
#'    Holds the structural zero probability of the ZIP distribution for each 
#'    observed count. These parameters are typically estimated from past data 
#'    using e.g. ZIP regression.
#' @param population Optional. A matrix or vector of populations for each 
#'    location. Only needed if \code{baselines} and \code{probs} are to be 
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
#' @param rel_tol A positive scalar. If the relative change in the incomplete
#'    information likelihood is less than this value, then the EM algorithm is
#'    deemed to have converged.
#' @return A list which, in addition to the information about the type of scan
#'    statistic, has the following components:
#'    \describe{
#'      \item{MLC}{A list containing the number of the zone of the most likely
#'            cluster (MLC), the locations in that zone, the duration of the 
#'            MLC, the calculated score, the relative risk, and the number of 
#'            iterations until convergence for the EM algorithm. In order, the 
#'            elements of this list are named  \code{zone_number, locations, 
#'            duration, score, relative_risk, n_iter}.}
#'      \item{observed}{A data frame containing, for each combination of zone 
#'            and duration investigated, the zone number, duration, score, 
#'            relative risk, number of EM iterations. The table is sorted by 
#'            score with the top-scoring location on top. If 
#'            \code{max_only = TRUE}, only contains a single row corresponding 
#'            to the MLC.}
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
#' @details For the expectation-based zero-inflated Poisson scan statistic
#'    (Allévius & Höhle 2017), the null hypothesis of no anomaly holds that 
#'    the count observed at each location \eqn{i} and duration \eqn{t} (the 
#'    number of time periods before present) has a zero-inflated Poisson 
#'    distribution with expected value parameter \eqn{\mu_{it}} and structural 
#'    zero probability \eqn{p_{it}}:
#'    \deqn{
#'      H_0 : Y_{it} \sim \textrm{ZIP}(\mu_{it}, p_{it}).
#'    }
#'    This holds for all locations \eqn{i = 1, \ldots, m} and all durations
#'    \eqn{t = 1, \ldots,T}, with \eqn{T} being the maximum duration considered.
#'    Under the alternative hypothesis, there is a space-time window \eqn{W}
#'    consisting of a spatial zone \eqn{Z \subset \{1, \ldots, m\}} and a time
#'    window \eqn{D \subseteq \{1, \ldots, T\}} such that the counts in that
#'    window have their Poisson expected value parameters inflated by a factor
#'    \eqn{q_W > 1} compared to the null hypothesis:
#'    \deqn{
#'    H_1 : Y_{it} \sim \textrm{ZIP}(q_W \mu_{it}, p_{it}), ~~(i,t) \in W.
#'    }
#'    For locations and durations outside of this window, counts are assumed to
#'    be distributed as under the null hypothesis. The sets \eqn{Z} considered
#'    are those specified in the argument \code{zones}, while the maximum
#'    duration \eqn{T} is taken as the maximum value in the column
#'    \code{duration} of the input \code{table}.
#'
#'    For each space-time window \eqn{W} considered, (the log of) a likelihood
#'    ratio is computed using the distributions under the alternative and null
#'    hypotheses, and the expectation-based Poisson scan statistic is calculated
#'    as the maximum of these quantities over all space-time windows. The
#'    expectation-maximization (EM) algorithm is used to obtain maximum
#'    likelihood estimates.
#' @references
#'    Allévius, B. and Höhle, M, \emph{An expectation-based space-time scan 
#'    statistic for ZIP-distributed data}, (Technical report),
#'    \href{https://goo.gl/yYJ42A}{Link to PDF}.
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
#' probs <- matrix(runif(n_total) / 4, max_duration, n_locs)
#' counts <- gamlss.dist::rZIP(n_total, baselines, probs)
#'
#' # Inject outbreak/event/anomaly
#' ob_dur <- 3
#' ob_cols <- zones[[10]]
#' ob_rows <- max_duration + 1 - seq_len(ob_dur)
#' counts[ob_rows, ob_cols] <- gamlss.dist::rZIP(
#'   ob_dur * length(ob_cols), 2 * baselines[ob_rows, ob_cols], 
#'   probs[ob_rows, ob_cols])
#' res <- scan_eb_zip(counts = counts,
#'                    zones = zones,
#'                    baselines = baselines,
#'                    probs = probs,
#'                    n_mcsim = 99,
#'                    max_only = FALSE,
#'                    rel_tol = 1e-3)
#' }
scan_eb_zip <- function(counts,
                        zones,
                        baselines = NULL,
                        probs = NULL,
                        population = NULL,
                        n_mcsim = 0,
                        max_only = FALSE,
                        rel_tol = 1e-3) {
  if (is.data.frame(counts)) {
    # Validate input -----------------------------------------------------------
    if (any(c("time", "location", "count") %notin% names(counts))) {
      stop("Data frame counts must have columns time, location, count")
    }
    counts %<>% arrange(location, -time)
    # Create matrices ----------------------------------------------------------
    if ("baseline" %in% names(counts)) {
      baselines <- df_to_matrix(counts, "time", "location", "baseline")
    }
    if ("prob" %in% names(counts)) {
      probs <- df_to_matrix(counts, "time", "location", "prob")
    }
    if ("population" %in% names(counts)) {
      population <- df_to_matrix(counts, "time", "location", "population")
    }
    counts <- df_to_matrix(counts, "time", "location", "count")
  }
  
  # Validate input -------------------------------------------------------------
  if (any(as.vector(counts) != as.integer(counts))) {
    stop("counts must be integer")
  }
  if (any(baselines <= 0)) stop("baselines must be positive")
  if (any(probs <= 0)) stop("probs must be positive")
  
  # Estimate baselines and probs if not supplied -------------------------------
  if (is.null(baselines) & is.null(population)) {
    stop("baselines or population matrices must be supplied")
  }
  if (is.null(baselines) || is.null(probs)) {
    warning("baselines and/or probs not supplied. ", 
            "Estimating ZIP parameters in heuristic fashion.")
    pars <- estimate_zip_params(counts, population)
    baselines <- pars$baselines
    probs <- pars$probs
  } 
  
  # Reshape into matrices ------------------------------------------------------
  if (is.vector(counts)) {
    counts <- matrix(counts, nrow = 1)
  }
  if (!is.null(baselines) && is.vector(baselines)) {
    baselines <- matrix(baselines, nrow = 1)
  }
  if (!is.null(probs) && is.vector(probs)) {
    probs <- matrix(probs, nrow = 1)
  }
  
  # Reverse time order: most recent first --------------------------------------
  counts <- flipud(counts)
  baselines <- flipud(baselines)
  probs <- flipud(probs)
  

  # Prepare zone arguments for C++ ---------------------------------------------
  args <- list(counts = counts, 
               baselines = baselines,
               probs = probs,
               zones = unlist(zones) - 1, 
               zone_lengths = unlist(lapply(zones, length)),
               rel_tol = rel_tol, 
               store_everything = !max_only,
               num_mcsim = n_mcsim)

  # Run analysis on observed counts --------------------------------------------
  scan <- run_scan(scan_eb_zip_cpp, args)

  MLC_row <- scan$observed[1, ]
  
  MLC_out <- list(zone_number = MLC_row$zone,
                  locations = zones[[MLC_row$zone]],
                  duration = MLC_row$duration,
                  score = MLC_row$score,
                  relative_risk = MLC_row$relrisk,
                  n_iter = MLC_row$n_iter)
  
  structure(
    c(list(# General
      distribution = "zero-inflated Poisson",
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
