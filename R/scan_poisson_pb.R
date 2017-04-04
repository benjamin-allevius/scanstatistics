# 
# #' @param counts
# #' @param population
# #' @keywords internal
# estimate_pop_baselines <- function(counts, population) {
#   # If analysis is purely spatial:
#   if (is.vector(counts)) {
#     if (!is.vector(population)) {
#       stop("If counts is a vector, population should be too.")
#     }
#     return(population * sum(counts) / sum(population))
#   }
#   
#   # If population is a vector, re-format to a matrix
#   if (is.vector(population)) {
#     population <- matrix(population, 
#                          nrow = nrow(counts),
#                          ncol = length(population),
#                          byrow = TRUE)
#   }
#   # Check that the same number of locations is implied by both arguments
#   if (ncol(counts) != ncol(population)) {
#     stop("The number of locations implied must be the same in the counts and ",
#          "the population arguments.")
#   }
#   # Else return baselines as fractions of the total count
#   timewise_pop <- rowSums(population)
#   baselines <- matrix(0, nrow(counts), ncol(counts))
#   for (i in seq_len(nrow(counts))) {
#     baselines[i, ] <- 
#       population[i, ] * (sum(counts) / (timewise_pop[i] * nrow(counts)))
#   }
#   return(baselines)
# }
# 
# #' Calculates the population-based Poisson LLR statistic for a given cluster.
# #' @param in_count A scalar; the observed count inside the cluster.
# #' @param in_base A scalar; the baseline (Poisson expected value) for the 
# #'    cluster.
# #' @param total_count A scalar; the total count over all locations and time 
# #'    points.
# #' @return A scalar.
# #' @keywords internal
# pb_poisson_statistic <- function(in_count, in_base, total_count) {
#   relrisk <- in_count / in_base
#   if (in_count <= in_base) {
#     return(c(score = 0, relative_risk = 1))
#   }
#   out_count <- total_count - in_count
#   out_base <- total_count - in_base
#   
#   term1 <- ifelse(in_count == 0, 0, in_count * (log(in_count / in_base)))
#   term2 <- ifelse(out_count == 0, 0, out_count * (log(out_count / out_base)))
#   return(c(score = term1 + term2, relative_risk = relrisk))
# }
# 
# calc_poisson_pb <- function(counts, baselines, zones, total_count = NULL) {
#   if (is.null(total_count)) {
#     total_count <- sum(counts)
#   }
#   n_dur <- nrow(counts)
#   
#   cs_counts <- apply(counts, 2, cumsum)
#   cs_base <- apply(baselines, 2, cumsum)
#   
#   result <- data.frame(zone = rep(seq_along(zones), n_dur),
#                        duration = rep(seq_len(n_dur), each = length(zones)),
#                        score = 0, 
#                        relative_risk = 0)
#   n_zones <- length(zones)
#   i <- seq_len(n_zones)
#   for (t in seq_len(n_dur)) {
#     result[i, c("score", "relative_risk")] <- 
#       t(sapply(zones, 
#                function(z) pb_poisson_statistic(cs_counts[1:t, z],
#                                                 cs_base[1:t, z],
#                                                 total_count),
#                simplify = "matrix"))
#   }
# }
# 
# #' @param counts A matrix: rows indicate time (top row most recent time point),
# #'    columns indicate locations (numbered 1 to \code{nrow(counts)}).
# #' @param population Either a matrix with the same dimensions as \code{counts},
# #'    or a vector of the same length as the number of columns in \code{counts}
# #'    (i.e. the number of locations).
# #' @param zones
# scan_poisson_pb <- function(counts, population, zones) {
#   
#   # If analysis is purely spatial:
#   if (is.vector(counts)) {
#     if (!is.vector(population)) {
#       stop("If counts is a vector, population should be too.")
#     }
#     counts <- matrix(counts, nrow = 1)
#     population <- matrix(population, nrow = 1)
#   }
#   
#   baselines <- estimate_pop_baselines(counts, population)
#   total_count <- sum(counts)
#   
#   
#   
#   
#   n_dur <- nrow(counts)
#   
#   
#   
#   
#   
#   
# }