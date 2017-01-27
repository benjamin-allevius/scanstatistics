
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
#'      \item{...}{Optional. Futher matrices with parameters}
#' }
#' @param score_fun A function taking matrix arguments, all of the
#'    same dimension, and returning a matrix or vector of that dimension. 
#' @param priority_fun A function taking matrix arguments, all of the
#'    same dimension, and returning a matrix or vector of that dimension. 
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
#' @inheritParams subset_aggregation
#' @param d An integer: \code{d=2} means sums are taken over locations (the 
#'    second dimension) and \code{d=3} means sums are taken over data streams.
#'    Other values of \code{d} should not be used.
#' @references 
#'    Neill, Daniel B., Edward McFowland, and Huanian Zheng (2013). \emph{Fast 
#'    subset scan for multivariate event detection}. Statistics in Medicine 
#'    32 (13), pp. 2185-2208.
#' @keywords internal
subset_aggregation <- function(args,
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
  subsets <- lapply(sets::set_power(seq_len(dim(counts)[d])), unlist)[-1]
  
  
  # Extract and sum over each data stream subset for each array in the list 
  # args. Then do Subset Aggregation with these sums.
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
  top_scoring
}
