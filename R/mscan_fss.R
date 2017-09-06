

#' Compute the Fast Subset Scan statistic for multivariate space-time data.
#' 
#' Compute the most likely cluster (MLC) using one of the Fast Subset Scan 
#' methods proposed by Neill et al. (2013). 
#' 
#' The data suitable for this function should consists of multiple variables 
#' ("data streams") observed over time at several locations, collected into an
#' array. The goal is to identify a subset (cluster) of data streams, locations, 
#' and time periods that has higher observed counts than expected. The method
#' will only detect clusters that are active, in the sense that they stretch 
#' from the most recent time period to some number of time periods back. The
#' counts can either be discrete or continuous; choose the \code{distribution}
#' parameter to suit your data.
#' 
#' @section Method: Subset Aggregation
#' Briefly, this method supposes the relative risk is constant and the same over
#' all data streams, locations, and time periods. Three versions of this method
#' exist, available through the parameter \code{algorithm}:
#' \describe{
#'   \item{fast}{Fast randomized optimization over both subsets of locations and 
#'               subsets of data streams.}
#'   \item{naive_streams}{Fast optimization over subsets of locations and naive 
#'                        optimization over subsets of streams. Can be used if 
#'                        the number of data streams is small. Denoted "FN" in
#'                        the paper by Neill et al. (2013).}
#'   \item{naive_locations}{Fast optimization over subsets of streams and naive 
#'                          optimization over subsets of locations. Can be used 
#'                          if the number of locations or spatial zones (groups
#'                          of locations considered jointly) is small. Denoted 
#'                          "NF" in the paper by Neill et al. (2013).}
#' }
#' Note: algorithm not quite as in Neill et al. (2013) since the randomly chosen 
#' subset of streams is the same for all time windows.
#' @section Method: Score Aggregation
#' Briefly, this method supposes that the relative risk is constant and the same
#' over all locations and time periods, but differ between data streams. Two 
#' versions of this method exist, available through the parameter 
#' \code{algorithm}:
#' \describe{
#'   \item{fast}{Fast randomized optimization over both subsets of locations and 
#'               subsets of data streams.}
#'   \item{naive_locations}{Fast optimization over subsets of streams and naive 
#'                          optimization over subsets of locations. Can be used 
#'                          if the number of locations or spatial zones (groups
#'                          of locations considered jointly) is small. Denoted 
#'                          "NK" in the paper by Neill et al. (2013).}
#' }
#' Note: this method is called "Kulldorff's method" in Neill et al. (2013).
#' @param counts An array of counts (integer or numeric). First
#'    dimension is time, ordered from most recent to most distant. Second 
#'    dimension indicates locations, which will be enumerated from 1 and up. 
#'    Third dimension indicates data streams, which will be enumerated from 1 
#'    and up.
#' @param distribution A string; one of "poisson", "gaussian", "exponential".
#' @param method A string; one of "subset" and "score". See explanation below.
#' @param algorithm A string; one of "fast", "naive_streams", "naive_locations".
#'    See explanation below.
#' @param parameters An optional list of parameters suitable for the 
#'    distribution chosen. Possible named elements are:
#'    \describe{
#'      \item{baselines}{An array of the same dimensions as \code{counts}. 
#'      Should hold the expected value of the count for each location, time
#'      point and data stream.}
#'      \item{variances}{An array of the same dimensions as \code{counts}.
#'      Should hold the variance of the count for each location, time
#'      point and data stream. Suitable for the gaussian distribution.}
#'    }
#' @param population An optional array, matrix or vector of populations.
#'    If an array, be of same dimensions as \code{counts}. If a matrix, should
#'    have as many rows as there are data streams and as many columns as there
#'    are locations. If a vector, should have the same length as the number of 
#'    locations. 
#' @param knn_matrix An optional integer matrix in which each row corresponds to
#'    a location. Each row starts with the index of the location (i.e. row 
#'    \eqn{i} has the integer \eqn{i} as its first element). Following that, the
#'    \code{(ncol(knn_matrix) - 1)} nearest neighbors of location \eqn{i} are
#'    listed in increasing order of distance on the same row. If this argument 
#'    is included, the search for the MLC are only done in these kNN subsets of
#'    locations.
#' @param zones An optional list of integer vectors. If included, the search for
#'    MLC will only be made in these subsets of locations.
#' @param ... Optional arguments, which are:
#'    \describe{
#'      \item{R}{The number of random restarts for the "fast" algorithms.}
#'      \item{rel_tol}{The relative tolerance criterion, used to determine 
#'      convergence for the "fast" algorithms. If the current score divided by 
#'      the previous score, minus one, is less than this number then the 
#'      algorithm is deemed to have converged.}
#'    }
#' @return A list containing the most likely cluster (MLC), having the following 
#'    elements:
#'    \describe{
#'      \item{score}{A scalar; the score of the MLC.}
#'      \item{duration}{An integer; the duration of the MLC, i.e. how many time 
#'                      periods from the present into the past the MLC 
#'                      stretches.}
#'      \item{locations}{An integer vector; the locations contained in the MLC.}
#'      \item{streams}{An integer vector; the data streams contained in the 
#'                     MLC.}
#'      \item{random_restarts}{FF only. The number of random restarts 
#'                             performed.}
#'      \item{iter_to_conv}{FF only. The number of iterations it took to reach 
#'                          convergence for each random restart.}
#'    }
#' @details 
#' @references 
#'    Neill, Daniel B., Edward McFowland, and Huanian Zheng (2013). \emph{Fast 
#'    subset scan for multivariate event detection}. Statistics in Medicine 
#'    32 (13), pp. 2185-2208.
#' @keywords internal
#' @examples 
#' \dontrun{
#' # Set simulation parameters (small)
#' set.seed(1)
#' n_loc <- 20
#' n_dur <- 10
#' n_streams <- 2
#' n_tot <- n_loc * n_dur * n_streams
#' 
#' # Create locations and kNN matrix
#' geo <- data.frame(x = rnorm(n_loc), y = rnorm(n_loc))
#' knn_mat <- coords_to_knn(geo, k = 10)
#' 
#' # Generate baselines and possibly other distribution parameters
#' baselines <- rexp(n_tot, 1/5) + rexp(n_tot, 1/5)
#' sigma2s <- rexp(n_tot)
#' 
#' # Generate counts
#' counts <- rpois(n_tot, baselines)
#' 
#' # Reshape into arrays
#' counts <- array(counts, c(n_dur, n_loc, n_streams))
#' baselines <- array(baselines, c(n_dur, n_loc, n_streams))
#' sigma2s <- array(sigma2s, c(n_dur, n_loc, n_streams))
#' 
#' # Inject an outbreak/event
#' ob_loc <- 1:floor(n_loc / 4)
#' ob_dur <- 1:floor(n_dur / 4)
#' ob_streams <- 1:floor(n_streams / 2)
#' counts[ob_dur, ob_loc, ob_streams] <- 4 * counts[ob_dur, ob_loc, ob_streams]
#' 
#' # Run the Subset Aggregation FN algorithm
#' FN_res <- mscan_fss(
#'   counts = counts,
#'   distribution = "poisson",
#'   algorithm = "naive_streams"
#'   parameters = list(baselines = baselines))
#'   
#' # Run the FF algorithm (few random restarts)
#' FF_res <- mscan_fss(
#'   counts = counts,
#'   distribution = "gaussian",
#'   algorithm = "fast"
#'   parameters = list(baselines = baselines, variances = variances),
#'   knn_matrix = knn_mat,
#'   R = 10)
#' }
mscan_fss <- function(counts,
                      distribution = c("poisson", "gaussian", "exponential"),
                      method = c("subset", "score"),
                      algorithm = c("fast", "naive_streams", "naive_locations"),
                      parameters = NULL,
                      population = NULL,
                      knn_matrix = NULL,
                      zones = NULL,
                      ...) {
  args <- list(...)
  args$counts <- counts
  
  distribution <- distribution[1]
  method <- method[1]
  algorithm <- algorithm[1]
  
  # Define score and priority functions matching the distribution
  if (distribution == "poisson") {
    priority_fun <- poisson_priority
    score_fun <- poisson_score
  } else if (distribution == "gaussian") {
    priority_fun <- gaussian_priority
    score_fun <- gaussian_score
  } else if (distribution == "exponential") {
    priority_fun <- exponential_priority
    score_fun <- exponential_score
  } else if ( !methods::hasArg("score_fun") || 
              !methods::hasArg("priority_fun") ) {
    stop("Either supply your own score and priority functions, or make sure",
         "that the argument distribution is one of 'poisson', 'gaussian'", 
         "'exponential'.")
  }
  
  # Define score and priority functions if they were supplied
  if (methods::hasArg("score_fun")) {
    score_fun <- args$score_fun
  }
  if (methods::hasArg("priority_fun")) {
    priority_fun <- args$priority_fun
  }
  
  # Estimate baselines and other parameters from data if not supplied
  if (is.null(parameters)) {
    args$baselines <- estimate_baselines(counts)
    
    if (distribution == "gaussian") {
      args$variances <- estimate_variances(counts)
    }
  }
  
  if (!(method %in% c("subset", "score"))) {
    stop('Method must be either "subset" or "score".')
  }
  
  # Subset aggregation ---------------------------------------------------------
  if (method == "subset") {
    if (algorithm == "fast") {
      return(subset_aggregation_FF(args,
                                   score_fun,
                                   priority_fun,
                                   ...))
    } else if (!(algorithm %in% c("naive_locations", "naive_streams"))) {
      stop('algorithm must be one of "fast", "naive_locations", ',
           '"naive_streams" if using the Subset Aggregation method.')
    } else {
      algo <- ifelse(algorithm == "naive_locations", "NF", "FN")
      return(subset_aggregation_FN_NF(args,
                                      score_fun,
                                      priority_fun,
                                      algo))
    }
  # Score aggregation ----------------------------------------------------------
  } else if (method == "score") {
    if (algorithm == "fast") {
      # return(score_aggregation_FF())
    } else if (algorithm != "naive_locations") {
      stop('algorithm must be one of "fast" and "naive_locations"')
    } else {
      # return(score_aggregation_NF())
    }
  } else {
    stop('method must be one of "subset" and "score".')
  }
}
