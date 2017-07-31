context("EB negbin statistic tests")

# Helper functions -----------------------------------------------------------
sc_hot <- function(y, m, w) {
  num <- sum((y - m) / w)
  den <- sum(m / w)
  return(num / den)
}

sc_inc <- function(y, m, w, d) {
  num <- vapply(1:ncol(y),
                function(i) sum( (y[ , i] - m[ , i]) * (d:1) / w[ , i]),
                numeric(1))
  den <- vapply(1:ncol(y),
                function(i) sum( m[ , i] * (d:1)^2 / w[ , i]),
                numeric(1))
  return(sum(num) / sum(den))
}

test_that("scan_eb_negbin: 1 timepoint, hotspot", {

  # Single timepoint -----------------------------------------------------------
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 2), nrow = 1),
    overdisp = matrix(c(1.5, 1.5), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))

  actual1 <- scan_eb_negbin_cpp(in1$counts,
                                in1$baselines,
                                in1$overdisp,
                                 in1$zones_flat - 1,
                                 in1$zone_lengths,
                                 store_everything = TRUE,
                                num_mcsim = 0,
                                score_hotspot = TRUE)$observed # hotspot
  actual1b <- scan_eb_negbin_cpp(in1$counts,
                                in1$baselines,
                                in1$overdisp,
                                in1$zones_flat - 1,
                                in1$zone_lengths,
                                store_everything = FALSE,
                                num_mcsim = 0,
                                score_hotspot = TRUE)$observed # hotspot

  expected1_score <- c(
    # Duration = 1
    sc_hot(in1$counts[1, 1],   in1$baselines[1, 1],   in1$overdisp[1, 1]),
    sc_hot(in1$counts[1, 2],   in1$baselines[1, 2],   in1$overdisp[1, 2]),
    sc_hot(in1$counts[1, 1:2], in1$baselines[1, 1:2], in1$overdisp[1, 1:2]))
  expect_equal(actual1$score, expected1_score)
  expect_equal(c(actual1[which.max(actual1$score), ]), c(actual1b))
})


test_that("scan_eb_negbin: 3 timepoints, hotspot", {

  # 3 timepoints ---------------------------------------------------------------
  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    baselines = matrix(c(0.5, 2,
                         0.5, 2,
                         0.5, 2), nrow = 3, byrow = TRUE),
    overdisp = matrix(c(1.5, 2,
                        1.5, 2,
                        1.5, 2), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))

  actual2 <- scan_eb_negbin_cpp(in2$counts,
                                in2$baselines,
                                in2$overdisp,
                                in2$zones_flat - 1,
                                in2$zone_lengths,
                                store_everything = TRUE,
                                num_mcsim = 0,
                                score_hotspot = TRUE)$observed # hotspot

  actual2b <- scan_eb_negbin_cpp(in2$counts,
                                in2$baselines,
                                in2$overdisp,
                                in2$zones_flat - 1,
                                in2$zone_lengths,
                                store_everything = FALSE,
                                num_mcsim = 0,
                                score_hotspot = TRUE)$observed # hotspot

  expected2_score <- c(
    # Duration = 1
    sc_hot(in2$counts[1, 1],   in2$baselines[1, 1],   in2$overdisp[1, 1]),
    sc_hot(in2$counts[1, 2],   in2$baselines[1, 2],   in2$overdisp[1, 2]),
    sc_hot(in2$counts[1, 1:2], in2$baselines[1, 1:2], in2$overdisp[1, 1:2]),
    # Duration = 2
    sc_hot(in2$counts[1:2, 1],   in2$baselines[1:2, 1],   in2$overdisp[1:2, 1]),
    sc_hot(in2$counts[1:2, 2],   in2$baselines[1:2, 2],   in2$overdisp[1:2, 2]),
    sc_hot(in2$counts[1:2, 1:2], in2$baselines[1:2, 1:2], in2$overdisp[1:2, 1:2]),
    # Duration = 3
    sc_hot(in2$counts[1:3, 1],   in2$baselines[1:3, 1],   in2$overdisp[1:3, 1]),
    sc_hot(in2$counts[1:3, 2],   in2$baselines[1:3, 2],   in2$overdisp[1:3, 2]),
    sc_hot(in2$counts[1:3, 1:2], in2$baselines[1:3, 1:2], in2$overdisp[1:3, 1:2]))

  expect_equal(actual2$score, expected2_score)
  expect_equal(c(actual2[which.max(actual2$score), ]), c(actual2b))
})


test_that("scan_eb_negbin: 1 timepoint, emerging", {

  # Single timepoint -----------------------------------------------------------
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 2), nrow = 1),
    overdisp = matrix(c(1.5, 1.5), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))

  actual1 <- scan_eb_negbin_cpp(in1$counts,
                                in1$baselines,
                                in1$overdisp,
                                in1$zones_flat - 1,
                                in1$zone_lengths,
                                store_everything = TRUE,
                                num_mcsim = 0,
                                score_hotspot = FALSE)$observed # emerging
  actual1b <- scan_eb_negbin_cpp(in1$counts,
                                 in1$baselines,
                                 in1$overdisp,
                                 in1$zones_flat - 1,
                                 in1$zone_lengths,
                                 store_everything = FALSE,
                                 num_mcsim = 0,
                                 score_hotspot = FALSE)$observed # emerging

  expected1_score <- c(
    # Duration = 1
    sc_inc(in1$counts[1, 1, drop=F],   in1$baselines[1, 1, drop=F],   in1$overdisp[1, 1, drop=F], 1),
    sc_inc(in1$counts[1, 2, drop=F],   in1$baselines[1, 2, drop=F],   in1$overdisp[1, 2, drop=F], 1),
    sc_inc(in1$counts[1, 1:2, drop=F], in1$baselines[1, 1:2, drop=F], in1$overdisp[1, 1:2, drop=F], 1))

  expect_equal(actual1$score, expected1_score)
  expect_equal(c(actual1[which.max(actual1$score), ]), c(actual1b))
})


test_that("scan_eb_negbin: 3 timepoints, emerging", {

  # 3 timepoints ---------------------------------------------------------------
  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    baselines = matrix(c(0.5, 2,
                         0.5, 2,
                         0.5, 2), nrow = 3, byrow = TRUE),
    overdisp = matrix(c(1.5, 2,
                        1.5, 2,
                        1.5, 2), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))

  actual2 <- scan_eb_negbin_cpp(in2$counts,
                                in2$baselines,
                                in2$overdisp,
                                in2$zones_flat - 1,
                                in2$zone_lengths,
                                store_everything = TRUE,
                                num_mcsim = 0,
                                score_hotspot = FALSE)$observed # hotspot

  actual2b <- scan_eb_negbin_cpp(in2$counts,
                                 in2$baselines,
                                 in2$overdisp,
                                 in2$zones_flat - 1,
                                 in2$zone_lengths,
                                 store_everything = FALSE,
                                 num_mcsim = 0,
                                 score_hotspot = FALSE)$observed # hotspot

  expected2_score <- c(
    # Duration = 1
    sc_inc(in2$counts[1, 1, drop=F],   in2$baselines[1, 1, drop=F],   in2$overdisp[1, 1, drop=F], 1),
    sc_inc(in2$counts[1, 2, drop=F],   in2$baselines[1, 2, drop=F],   in2$overdisp[1, 2, drop=F], 1),
    sc_inc(in2$counts[1, 1:2, drop=F], in2$baselines[1, 1:2, drop=F], in2$overdisp[1, 1:2, drop=F], 1),
    # Duration = 2
    sc_inc(in2$counts[1:2, 1, drop=F],   in2$baselines[1:2, 1, drop=F],   in2$overdisp[1:2, 1, drop=F], 2),
    sc_inc(in2$counts[1:2, 2, drop=F],   in2$baselines[1:2, 2, drop=F],   in2$overdisp[1:2, 2, drop=F], 2),
    sc_inc(in2$counts[1:2, 1:2, drop=F], in2$baselines[1:2, 1:2, drop=F], in2$overdisp[1:2, 1:2], 2),
    # Duration = 3
    sc_inc(in2$counts[1:3, 1, drop=F],   in2$baselines[1:3, 1, drop=F],   in2$overdisp[1:3, 1, drop=F], 3),
    sc_inc(in2$counts[1:3, 2, drop=F],   in2$baselines[1:3, 2, drop=F],   in2$overdisp[1:3, 2, drop=F], 3),
    sc_inc(in2$counts[1:3, 1:2, drop=F], in2$baselines[1:3, 1:2, drop=F], in2$overdisp[1:3, 1:2, drop=F], 3))

  expect_equal(actual2$score, expected2_score)
  expect_equal(c(actual2[which.max(actual2$score), ]), c(actual2b))
})
