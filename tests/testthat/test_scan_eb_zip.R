context("EB ZIP statistic tests")



test_that("scan_eb_zip_cpp", {

  # Helper functions
  zip_lpmf <- function(x, mu, p) {
    if (x == 0) return(log(p + (1-p) * exp(-mu)))
    return(log(1-p) + dpois(x, mu, log=TRUE))
  }

  est_eb_zip_zeroindic <- function(mu, p, q) p / (p + (1 - p) * exp(-q * mu))

  zip_loglihood <- function(y, mu, p, q) {
    loglihood = 0
    for (i in 1:length(y)) {
      loglihood <- loglihood + zip_lpmf(y[i], q * mu[i], p[i])
    }
    return(loglihood)
  }


  # Single timepoint
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 2), nrow = 1),
    probs = matrix(c(0.1, 0.2), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))

  actual1 <- scan_eb_zip_cpp(in1$counts,
                             in1$baselines,
                             in1$probs,
                             in1$zones_flat - 1,
                             in1$zone_lengths,
                             rel_tol = 1e-3,
                             store_everything = TRUE,
                             num_mcsim = 0)$observed
  actual1b <- scan_eb_zip_cpp(in1$counts,
                             in1$baselines,
                             in1$probs,
                             in1$zones_flat - 1,
                             in1$zone_lengths,
                             rel_tol = 1e-3,
                             store_everything = FALSE,
                             num_mcsim = 0)$observed

  expected1_score <- c(zip_lpmf(1, 1, 0.1) - zip_lpmf(1, 0.5, 0.1), 0, 0)

  expect_equal(actual1$score, expected1_score)
  expect_equal(actual1$relrisk, c(2, 1, 1))
  expect_equal(actual1[which.max(actual1$score), ], actual1b)

  # 3 timepoints
  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    baselines = matrix(c(0.5, 2,
                         0.5, 2,
                         0.5, 2), nrow = 3, byrow = TRUE),
    probs = matrix(c(0.1, 0.2,
                     0.1, 0.2,
                     0.1, 0.2), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))

  actual2 <- scan_eb_zip_cpp(in2$counts,
                             in2$baselines,
                             in2$probs,
                             in2$zones_flat - 1,
                             in2$zone_lengths,
                             rel_tol = 1e-3,
                             store_everything = TRUE,
                             num_mcsim = 0)$observed
  actual2b <- scan_eb_zip_cpp(in2$counts,
                              in2$baselines,
                              in2$probs,
                              in2$zones_flat - 1,
                              in2$zone_lengths,
                              rel_tol = 1e-3,
                              store_everything = FALSE,
                              num_mcsim = 0)$observed

  recursive_zic <- function(n, y, mu, p) {
    mu1 <- sum(mu[y != 0])
    mu0 <- mu[y == 0]
    p0 <- p[y == 0]
    q <- 1
    for (i in 1:n) {
      d <- numeric(length(mu0))
      for (j in 1:length(mu0)) {
        d[j] <- est_eb_zip_zeroindic(mu0[j], p0[j], q)
      }
      q <- sum(y) / (mu1 + sum((1-d) * mu0))
    }
    return(q)
  }

  expected2_relrisk <- c(
    2,
    1,
    1,
    3,
    1,
    recursive_zic(5, in2$counts[1:2, 1:2], in2$baselines[1:2, 1:2], in2$probs[1:2, 1:2]),
    recursive_zic(5, in2$counts[1:3, 1], in2$baselines[1:3, 1], in2$probs[1:3, 1]),
    recursive_zic(5, in2$counts[1:3, 2], in2$baselines[1:3, 2], in2$probs[1:3, 2]),
    recursive_zic(5, in2$counts, in2$baselines, in2$probs))

  expected2_score <- c(
    # Duration = 1
    zip_lpmf(1, expected2_relrisk[1] * 0.5, 0.1) - zip_lpmf(1, 0.5, 0.1),
     0, 0,
    # Duration = 2
     zip_lpmf(3, expected2_relrisk[4] * 1, 0.1) - zip_lpmf(3, 1, 0.1),
     0,
     zip_loglihood(c(1, 2, 0, 1), c(0.5, 0.5, 2, 2), c(0.1, 0.1, 0.2, 0.2),
                   expected2_relrisk[6]) -
       zip_loglihood(c(1, 2, 0, 1), c(0.5, 0.5, 2, 2), c(0.1, 0.1, 0.2, 0.2),
                     1),
    # Duration = 3
    # zone = 1
     zip_loglihood(c(1, 2, 0), rep(0.5, 3), rep(0.1, 3), expected2_relrisk[7]) -
       zip_loglihood(c(1, 2, 0), rep(0.5, 3), rep(0.1, 3), 1),
    # zone = 2
     zip_loglihood(c(0, 1, 20), rep(2, 3), rep(0.2, 3), expected2_relrisk[8]) -
       zip_loglihood(c(0, 1, 20), rep(2, 3), rep(0.2, 3), 1),
    zip_loglihood(in2$counts, in2$baselines, in2$probs, expected2_relrisk[9]) -
      zip_loglihood(in2$counts, in2$baselines, in2$probs, 1))

  expect_equal(actual2$score, expected2_score, tolerance = 1e-3)
  expect_equal(actual2$relrisk, expected2_relrisk, tolerance = 1e-3)
  expect_equal(c(actual2[which.max(actual2$score), ]), c(actual2b))
})
