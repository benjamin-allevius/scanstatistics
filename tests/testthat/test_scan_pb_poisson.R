# context("PB Poisson statistic tests")
# 
# test_that("scan_pb_poisson_cpp", {
#   # Single timepoint
#   in1 <- list(
#     counts = matrix(c(1, 0), nrow = 1),
#     baselines = matrix(c(0.5, 0.5), nrow = 1),
#     zones = list(1L, 2L, 1:2))
#   in1$N = sum(in1$counts)
#   in1$zones_flat =  unlist(in1$zones)
#   in1$zone_lengths = unlist(lapply(in1$zones, length))
#   
#   actual1 <- scan_pb_poisson_cpp(
#     in1$counts, in1$baselines, in1$zones_flat - 1, in1$zone_lengths, in1$N)
#   expected1_score <- c(log(1 / 0.5), 0, 0)
#   expect_equal(actual1$score, expected1_score)
#   expect_equal(actual1$relrisk_in, c(2, 0, 1))
#   expect_equal(actual1$relrisk_out, c(0, 2, 1))
#   
#   # 3 timepoints
#   in2 <- list(
#     counts = matrix(c(1, 0,
#                       2, 1,
#                       0, 20), nrow = 3, byrow = TRUE),
#     baselines = matrix(c(2, 1.0,
#                          4, 5,
#                          6, 6), nrow = 3, byrow = TRUE),
#     zones = list(1L, 2L, 1:2))
#   in2$zones_flat =  unlist(in2$zones)
#   in2$zone_lengths = unlist(lapply(in2$zones, length))
#   in2$N <- sum(in2$counts)
#   
#   actual2 <- scan_pb_poisson_cpp(in2$counts,
#                                  apply(in2$baselines, 2, cumsum),
#                                  in2$zones_flat - 1, in2$zone_lengths,
#                                  in2$N)
#   
#   expected2_relrisk_in <- c(1/2, 0, 1/3,
#                             1/2, 1/6, 4/12,
#                             3 / 12, 21 / 12, 1)
#   expected2_relrisk_out <- c(23 / 22, 24 / 23, 23 / 21,
#                              21 / 18, 23 / 18, 20 / 12,
#                              21 / 12, 3 / 12, 1)
#   
#   f <- function(c, b) {
#     term2 <- ifelse(in2$N > b, (in2$N - c) * log((in2$N - c) / (in2$N - b)), 0)
#     ifelse(c > b, c * log(c/b) + term2, 0)
#   }
#   expected2_score <- c(f(1, 2), f(0, 1), f(1, 3),
#                        f(3, 6), f(1, 6), f(4, 12),
#                        f(3, 12), f(21, 12), f(24, 24))
#   
#   expect_equal(actual2$score, expected2_score)
#   expect_equal(actual2$relrisk_in, expected2_relrisk_in)
#   expect_equal(actual2$relrisk_out, expected2_relrisk_out)
# })
