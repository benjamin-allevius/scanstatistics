context("EB Poisson statistic tests")

test_that("scan_eb_poisson", {

  # Helper functions -----------------------------------------------------------
  poisson_lpmf <- function(y, mu) -mu + y * log(mu)
  poisson_loglihood <- function(y, mu, q) {
    llh <- 0
    for (i in 1:length(y)) {
      llh <- llh + poisson_lpmf(y[i], q * mu[i])
    }
    return(llh)
  }

  # Single timepoint -----------------------------------------------------------
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 2), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))

  actual1 <- scan_eb_poisson_cpp(in1$counts,
                                 in1$baselines,
                                 in1$zones_flat - 1,
                                 in1$zone_lengths,
                                 store_everything = TRUE,
                                 num_mcsim = 0)$observed
  actual1b <- scan_eb_poisson_cpp(in1$counts,
                                 in1$baselines,
                                 in1$zones_flat - 1,
                                 in1$zone_lengths,
                                 store_everything = FALSE,
                                 num_mcsim = 0)$observed
  expected1_score <- c(poisson_lpmf(1, 1) - poisson_lpmf(1, 0.5), 0, 0)
  expect_equal(actual1$score, expected1_score)
  expect_equal(actual1$relrisk, c(2, 1, 1))
  expect_equal(c(actual1[which.max(actual1$score), ]), c(actual1b))

  # 3 timepoints ---------------------------------------------------------------
  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    baselines = matrix(c(0.5, 2,
                         0.5, 2,
                         0.5, 2), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))

  actual2 <- scan_eb_poisson_cpp(in2$counts,
                                 in2$baselines,
                                 in2$zones_flat - 1,
                                 in2$zone_lengths,
                                 store_everything = TRUE,
                                 num_mcsim = 0)$observed
  actual2b <- scan_eb_poisson_cpp(in2$counts,
                                 in2$baselines, 
                                 in2$zones_flat - 1,
                                 in2$zone_lengths,
                                 store_everything = FALSE,
                                 num_mcsim = 0)$observed

  expected2_relrisk <- c(# Duration = 1
                         sum(in2$counts[1, 1]) / sum(in2$baselines[1, 1]),
                         sum(in2$counts[1, 2]) / sum(in2$baselines[1, 2]),
                         sum(in2$counts[1, 1:2]) / sum(in2$baselines[1, 1:2]),
                         # Duration = 2
                         sum(in2$counts[1:2, 1]) / sum(in2$baselines[1:2, 1]),
                         sum(in2$counts[1:2, 2]) / sum(in2$baselines[1:2, 2]),
                         sum(in2$counts[1:2, 1:2]) / sum(in2$baselines[1:2, 1:2]),
                         # Duration = 3
                         sum(in2$counts[1:3, 1]) / sum(in2$baselines[1:3, 1]),
                         sum(in2$counts[1:3, 2]) / sum(in2$baselines[1:3, 2]),
                         sum(in2$counts[1:3, 1:2]) / sum(in2$baselines[1:3, 1:2]))
  expected2_relrisk <- pmax(expected2_relrisk, 1)
  relrisk_mat <- matrix(expected2_relrisk, nrow = 3, ncol = 3, byrow = TRUE)

  expected2_score <- c(
    # Duration = 1
    poisson_lpmf(1, expected2_relrisk[1] * 0.5) - poisson_lpmf(1, 0.5), 0, 0,
    # Duration = 2
    poisson_lpmf(3, expected2_relrisk[4] * 1) - poisson_lpmf(3, 1), 0, 0,
    # Duration = 3
    # zone = 1
     poisson_loglihood(c(1, 2, 0), rep(0.5, 3), expected2_relrisk[7]) -
       poisson_loglihood(c(1, 2, 0), rep(0.5, 3), 1),
    # zone = 2
     poisson_loglihood(c(0, 1, 20), rep(2, 3), expected2_relrisk[8]) -
       poisson_loglihood(c(0, 1, 20), rep(2, 3), 1),
    # zone = 3
    poisson_loglihood(as.vector(in2$counts), as.vector(in2$baselines),
                      expected2_relrisk[9]) -
      poisson_loglihood(as.vector(in2$counts), as.vector(in2$baselines), 1))

  expect_equal(actual2$relrisk, expected2_relrisk)
  expect_equal(actual2$score, expected2_score)
  expect_equal(c(actual2[which.max(actual2$score), ]), c(actual2b))
})
