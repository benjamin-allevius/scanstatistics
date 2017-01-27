
#' @param counts A matrix of counts. Rows indicate time, ordered from most 
#'    recent to most distant. Columns indicate e.g. locations or data streams, 
#'    enumerated from 1 and up.
#' @param baselines A matrix of expected counts. Dimensions are as for 
#'    \code{counts}.
#' @param penalties A matrix of penalty terms. Dimensions are as for 
#'    \code{counts}.
#' @param score_fun A function taking matrix arguments, all of the
#'    same dimension, and returning a matrix or vector of that dimension. 
#' @param priority_fun A function taking matrix arguments, all of the
#'    same dimension, and returning a matrix or vector of that dimension. 
#' @param ... Arguments passed to \code{score_fun} and/or \code{priority_fun}. 
#'    Must be matrices of the same dimension as \code{counts} and 
#'    \code{baselines}.
#' @return A list containing three elements:
#'    \describe{
#'      \item{score}{The highest score of all clusters.}
#'      \item{duration}{The duration of the score-maximizing cluster.}
#'      \item{subset}{An integer vector of the subset of e.g. locations or data
#'                    streams in the score-maximizing cluster.}
#' }
#' @keywords internal
algo1 <- function(counts, 
                  baselines, 
                  penalties = NULL,
                  score_fun = poisson_score,
                  priority_fun = poisson_priority,
                  ...) {
  
  # List of matrices with rows=time and columns=locations/data streams
  args <- list(counts = counts, baselines = baselines, ...)
  
  # Compute location/data stream priorities and sort them thereafter
  prios <- do.call(priority_fun, args)
  prio_indices <- prioritize_cols(prios)
  
  args$priority_indices <- prio_indices
  
  #
  scores <- do.call(score_fun, args)
  
  if (!is.null(penalties)) {
    scores <- scores + reorder_rows(apply(penalties, 2, cumsum), prio_indices)
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

FN_SA <- function(counts, 
                  baselines, 
                  score_fun = function(c, b) ifelse(c > b, c * log(c/b) + b - c, 0)) {
  # With n_streams data streams labelled 1,...,n_streams, form a list of the 
  # 2^(n_streams)-1 non-empty subsets of the integers {1,...,n_streams}
  stream_subsets <- lapply(sets::set_power(seq_len(dim(counts)[3])), unlist)[-1]
  
  # For each stream subset, extract the corresponding columns in the 
  # data/parameter arrays, and perform algorithm 1
  res <- lapply(stream_subsets, 
                function(x) algo1(counts[ , , x], baselines[ , , x], score_fun))
  
  # Extract scores and locate the top-scoring location-duration subset (MLC)
  scores <- vapply(res, function(lst) lst[["score"]], numeric(1))
  maximizer <- which.max(scores)
  
  # Add the stream subset to the MLC and return
  top_scoring <- res[[maximizer]]
  top_scoring$streams <- stream_subsets[[maximizer]]
  top_scoring
}
