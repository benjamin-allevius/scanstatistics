# Fast Subset Score

#' Compute the highest scoring subset for aggregated data.
#' 
#' Given data that is aggregated over either locations or data streams, this
#' function finds the highest scoring subset of the remaining dimensions of the
#' data. For example, if the data has been aggregated over data streams, the
#' highest scoring subset will consist of a window of time stretching from the
#' most recent time period to some time period in the past (i.e. a 
#' \emph{duration}) and a collection of locations.
#' @param args A list of matrices:
#'    \describe{
#'      \item{counts}{Required. A matrix of counts. Rows indicate time, ordered 
#'                    from most recent to most distant. Columns indicate e.g. 
#'                    locations or data streams, enumerated from 1 and up.}
#'      \item{baselines}{Required. A matrix of expected counts. Dimensions are 
#'                       as for \code{counts}.}
#'      \item{penalties}{Optional. A matrix of penalty terms. Dimensions are as
#'                       for \code{counts}.}
#'      \item{...}{Optional. More matrices with parameters}
#' }
#' @param score_fun A function taking matrix arguments, all of the
#'    same dimension, and returning a matrix or vector of that dimension. 
#'    Suitable alternatives are \code{\link{poisson_score}}, 
#'    \code{\link{gaussian_score}}.
#' @param priority_fun A function taking matrix arguments, all of the
#'    same dimension, and returning a matrix or vector of that dimension. 
#'    Suitable alternatives are \code{\link{poisson_priority}}, 
#'    \code{\link{gaussian_priority}}.
#' @return A list containing three elements:
#'    \describe{
#'      \item{score}{The highest score of all clusters.}
#'      \item{duration}{The duration of the score-maximizing cluster.}
#'      \item{subset}{An integer vector of the subset of e.g. locations or data
#'                    streams in the score-maximizing cluster.}
#' }
#' @details This function provides the main component of the \emph{FN} and 
#'    \emph{NF} algorithms described in Section 3.1 of Neill et al. (2013).
#' @references 
#'    Neill, Daniel B., Edward McFowland, and Huanian Zheng (2013). \emph{Fast 
#'    subset scan for multivariate event detection}. Statistics in Medicine 
#'    32 (13), pp. 2185-2208.
#' @keywords internal
score_priority_subset <- function(args,
                                  score_fun = poisson_score,
                                  priority_fun = poisson_priority) {
  
  # Compute location/data stream priorities and sort them thereafter
  prios <- do.call(priority_fun, args)
  prio_indices <- prioritize_cols(prios)
  
  args$priority_indices <- prio_indices
  
  #
  scores <- do.call(score_fun, args)
  
  if ("penalties" %in% names(args)) {
    scores <- scores + reorder_rows(apply(args$penalties, 2, cumsum), 
                                    prio_indices)
  }
  
  # Get the index of the row (duration) and column (location subset) that 
  # maximizes the score
  max_score <- max(scores)
  maxer <- which(scores == max_score, arr.ind = TRUE)
  duration <- unname(maxer[1, 1])
  subset <- prio_indices[duration, seq_len(maxer[1, 2])]
  
  list(duration = duration,
       subset = subset,
       score = max_score)
}

#' Compute the most likely cluster using the FN/NF Subset Aggregation algorithm.
#' 
#' Compute the most likely cluster (MLC) with the Subset Aggregation method by
#' Neill et al. (2013), either through fast optimization over subsets of 
#' locations and naive optimization over subsets of streams (FN), or through
#' naive optimization over subsets of locations and fast optimization over 
#' subsets of streams (NF).
#' @param args A list of arrays:
#'    \describe{
#'      \item{counts}{Required. An array of counts (integer or numeric). First
#'                    dimension is time, ordered from most recent to most 
#'                    distant. Second dimension indicates locations, which will 
#'                    be enumerated from 1 and up. Third dimension indicates 
#'                    data streams, which will be enumerated from 1 and up.}
#'      \item{baselines}{Required. A matrix of expected counts. Dimensions are 
#'                       as for \code{counts}.}
#'      \item{penalties}{Optional. A matrix of penalty terms. Dimensions are as
#'                       for \code{counts}.}
#'      \item{...}{Optional. More matrices with distribution parameters.
#'                 Dimensions are as for \code{counts}.}
#' }
#' @inheritParams score_priority_subset
#' @param algorithm Either "FN" or "NF":
#'    \describe{
#'      \item{FN}{Fast optimization over subsets of locations and naive 
#'                optimization over subsets of streams. Can be used if the 
#'                number of data streams is small.}
#'      \item{NF}{Fast optimization over subsets of streams and naive 
#'                optimization over subsets of locations. Can be used if the 
#'                number of locations is small.}
#'    }
#' @return A list with 4 elements:
#'    \describe{
#'      \item{score}{A scalar; the score of the MLC.}
#'      \item{duration}{An integer; the duration of the MLC, i.e. how many time 
#'                      periods from the present into the past the MLC 
#'                      stretches.}
#'      \item{locations}{An integer vector; the locations contained in the MLC.}
#'      \item{streams}{An integer vector; the data streams contained in the 
#'                     MLC.}
#'   }
#' @references 
#'    Neill, Daniel B., Edward McFowland, and Huanian Zheng (2013). \emph{Fast 
#'    subset scan for multivariate event detection}. Statistics in Medicine 
#'    32 (13), pp. 2185-2208.
#' @keywords internal
subset_aggregation_FN_NF <- function(args,
                                     score_fun = poisson_score,
                                     priority_fun = poisson_priority,
                                     algorithm = "FN") {
  if (algorithm == "FN") {
    d <- 3
    fast_name <- "locations"
    naive_name <- "streams"
  } else if (algorithm == "NF") {
    d <- 2
    fast_name <- "streams"
    naive_name <- "locations"
  } else {
    stop("algorithm must be FN or NF")
  }
  
  # Naive optimization: iterate over power set (with empty set removed)
  subsets <- lapply(sets::set_power(seq_len(dim(args$counts)[d])), unlist)[-1]
  
  # Extract and sum over each data stream/location subset for each array in the 
  # list args. Then do Subset Aggregation with these sums.
  res <- lapply(subsets, 
                function(x) score_priority_subset(sum_over_subset(args, x, d),
                                                  score_fun,
                                                  priority_fun))
  
  # Extract scores and locate the top-scoring location-duration subset (MLC)
  scores <- vapply(res, function(lst) lst[["score"]], numeric(1))
  maximizer <- which.max(scores)
  
  # Add the location/stream subset to the MLC and return
  top_scoring <- res[[maximizer]]
  top_scoring[[naive_name]] <- subsets[[maximizer]]
  names(top_scoring)[which(names(top_scoring) == "subset")] <- fast_name
  top_scoring[c(3, 1, 2, 4)] # Return in same order as subset_aggregation_FF
}

#' Fast Subset Aggregation over both locations and data streams.
#' 
#' Compute the most likely cluster (MLC) with the Subset Aggregation method by
#' Neill et al. (2013) through fast optimization over subsets of locations and 
#' subsets of streams.
#' @inheritParams score_priority_subset
#' @param R The number of random restarts.
#' @param rel_tol The relative tolerance criterion. If the current score divided
#'    by the previous score, minus one, is less than this number then the 
#'    algorithm is deemed to have converged.
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
#'      \item{random_restarts}{The number of random restarts performed.}
#'      \item{iter_to_conv}{The number of iterations it took to reach 
#'                          convergence for each random restart.}
#'    }
#' @details Note: algorithm not quite as in Neill et al. (2013) since the 
#'    randomly chosen subset of streams is the same for all time windows.
#' @references 
#'    Neill, Daniel B., Edward McFowland, and Huanian Zheng (2013). \emph{Fast 
#'    subset scan for multivariate event detection}. Statistics in Medicine 
#'    32 (13), pp. 2185-2208.
#' @importFrom stats rbinom runif
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
#' # Run the FN algorithm
#' FF_res <- subset_aggregation_FF(
#'   list(counts = counts, baselines = baselines),
#'   score_fun = poisson_score,
#'   priority_fun = poisson_priority,
#'   algorithm = "FN")
#' }
subset_aggregation_FF <- function(args,
                                  score_fun = poisson_score,
                                  priority_fun = poisson_priority,
                                  R = 50,
                                  rel_tol = 1e-2) {
  dims <- dim(args[[1]])
  n_streams <- dims[3]
  
  all_scores <- numeric(R)
  results <- vector("list", R) 
  
  n_iterations <- numeric(R)
  
  i <- 1
  while (i <= R) {
    
    # Pick a random subset of locations
    stream_subset <- seq_len(n_streams)[rbinom(n_streams, 1, runif(1)) == 1]
    if (length(stream_subset) == 0) {
      next
    }
    score_prev <- 1
    score <- (1 + 2 * rel_tol) * score_prev
    
    iterations <- 1
    
    while (abs((score - score_prev) / score_prev) > rel_tol) {
      
      score_prev <- score
      
      # Perform the algorithm at the core of the FN method and get the highest
      # scoring subset of locations, conditional on the set of streams
      FN_score <- score_priority_subset(sum_over_subset(args, stream_subset, 3),
                                        score_fun,
                                        priority_fun)
      loc_subset <- FN_score$subset
      
      # Given the conditionally optimal subset of locations, perform the 
      # algorithm at the core of the NF method and get the highest
      # scoring subset of streams
      NF_score <- score_priority_subset(sum_over_subset(args, loc_subset, 2),
                                        score_fun,
                                        priority_fun)
      stream_subset <- NF_score$subset
      
      score <- NF_score$score
      
      iterations <- iterations + 1
      
    }
    
    all_scores[i] <- score
    results[[i]] <- list(score = score,
                         duration = NF_score$duration,
                         locations = loc_subset,
                         streams = stream_subset,
                         random_restarts = R)
    n_iterations[i] <- iterations
    
    i <- i + 1
  }
  
  res <- results[[which.max(all_scores)]]
  res$iter_to_conv <- n_iterations
  res
}

#' Subset Aggregation over locations and data streams, naive or fast.
#' 
#' Compute the most likely cluster (MLC) using either of three versions of the
#' Subset Aggregation method by Neill et al. (2013). The methods are:
#' \describe{
#'   \item{FF}{Fast optimization over both subsets of locations and subsets 
#'             of data streams.}
#'   \item{FN}{Fast optimization over subsets of locations and naive 
#'             optimization over subsets of streams. Can be used if the 
#'             number of data streams is small.}
#'   \item{NF}{Fast optimization over subsets of streams and naive 
#'             optimization over subsets of locations. Can be used if the 
#'             number of locations is small.}
#' }
#' @inheritParams subset_aggregation_FN_NF
#' @param R The number of random restarts.
#' @param rel_tol The relative tolerance criterion. If the current score divided
#'    by the previous score, minus one, is less than this number then the 
#'    algorithm is deemed to have converged.
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
#' @details Note: algorithm not quite as in Neill et al. (2013) since the 
#'    randomly chosen subset of streams is the same for all time windows.
#' @references 
#'    Neill, Daniel B., Edward McFowland, and Huanian Zheng (2013). \emph{Fast 
#'    subset scan for multivariate event detection}. Statistics in Medicine 
#'    32 (13), pp. 2185-2208.
#' @importFrom stats rbinom runif
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
#' # Run the FN algorithm
#' FN_res <- subset_aggregation(
#'   list(counts = counts, baselines = baselines),
#'   score_fun = poisson_score,
#'   priority_fun = poisson_priority,
#'   algorithm = "FN")
#'   
#' # Run the FF algorithm (few random restarts)
#' FN_res <- subset_aggregation(
#'   list(counts = counts, baselines = baselines),
#'   score_fun = poisson_score,
#'   priority_fun = poisson_priority,
#'   algorithm = "FN",
#'   R = 10)
#' }
subset_aggregation <- function(args,
                               score_fun = poisson_score,
                               priority_fun = poisson_priority,
                               algorithm = "FF",
                               R = 50,
                               rel_tol = 1e-2) {
  if (algorithm == "FF") {
    return(subset_aggregation_FF(args,
                                 score_fun,
                                 priority_fun,
                                 R,
                                 rel_tol))
  } else if (!(algorithm %in% c("FN", "NF"))) {
    stop('algorithm must be one of "FF", "FN", "NF".')
  } else {
    return(subset_aggregation_FN_NF(args,
                                    score_fun,
                                    priority_fun,
                                    algorithm))
  }
}
