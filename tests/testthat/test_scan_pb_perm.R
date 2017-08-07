context("PB Poisson statistic tests")

test_that("scan_pb_perm_cpp", {

  # 3 timepoints
  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$baselines <- outer(rowSums(A), colSums(A)) / sum(A)
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))
  in2$N <- sum(in2$counts)

  actual2 <- scan_pb_perm_cpp(in2$counts,
                              in2$baselines,
                              in2$zones_flat - 1,
                              in2$zone_lengths,
                              store_everything = TRUE,
                              num_mcsim = 0)$observed
  actual2b <- scan_pb_poisson_cpp(in2$counts,
                                  in2$baselines,
                                 in2$zones_flat - 1,
                                 in2$zone_lengths,
                                 store_everything = FALSE,
                                 num_mcsim = 0)$observed

  expected2_relrisk_in <- c(1/2, 0, 1/3,
                            1/2, 1/6, 4/12,
                            3 / 12, 21 / 12, 1)
  expected2_relrisk_out <- c(23 / 22, 24 / 23, 23 / 21,
                             21 / 18, 23 / 18, 20 / 12,
                             21 / 12, 3 / 12, 1)

  f <- function(c, b) {
    term2 <- ifelse(in2$N > b, (in2$N - c) * log((in2$N - c) / (in2$N - b)), 0)
    ifelse(c > b, c * log(c/b) + term2, 0)
  }
  expected2_score <- c(f(1, 2), f(0, 1), f(1, 3),
                       f(3, 6), f(1, 6), f(4, 12),
                       f(3, 12), f(21, 12), f(24, 24))

  expect_equal(actual2$score, expected2_score)
  expect_equal(actual2$relrisk_in, expected2_relrisk_in)
  expect_equal(actual2$relrisk_out, expected2_relrisk_out)
  expect_equal(c(actual2[which.max(actual2$score), ]), c(actual2b))
})

test_that("scan_pb_poisson_cpp mcsim", {

  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    baselines = matrix(c(2, 1.0,
                         4, 5,
                         6, 6), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))
  in2$N <- sum(in2$counts)

  # Generate counts with rnbinom with R, run scan_pb_poisson_cpp
  set.seed(1)
  num_sims <- 2
  obs_res <- NULL
  for (i in 1:num_sims) {
    in2$counts <- matrix(rmultinom(1, in2$N, as.vector(in2$baselines)),
                       nrow(in2$baselines), ncol(in2$baselines))
    res <- scan_pb_poisson_cpp(in2$counts,
                               in2$baselines,
                               in2$zones_flat - 1,
                               in2$zone_lengths,
                               store_everything = FALSE,
                               num_mcsim = 0)$observed
    obs_res <- rbind(obs_res, res)
  }

  # Run equal number of mcsim
  set.seed(1)
  sim_res <- scan_pb_poisson_cpp(in2$counts,
                                in2$baselines,
                                in2$zones_flat - 1,
                                in2$zone_lengths,
                                store_everything = FALSE,
                                num_mcsim = 2)$simulated

  expect_equal(obs_res, sim_res)
})
