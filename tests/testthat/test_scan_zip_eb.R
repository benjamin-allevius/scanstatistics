context("EB ZIP statistic tests")


test_that("est_eb_zip_relrisk", {
  y1 <- c(1, 0); mu1 <- c(0.5, 2); p1 <- c(0.1, 0.2); d1 <- c(0, 0.65)
  expect_equal(est_eb_zip_relrisk(sum(y1), mu1, p1, d1), 1)
  
  y2 <- c(2, 2); mu2 <- c(1, 1); p2 <- c(0.1, 0.2); d2 <- c(0, 0)
  expect_equal(est_eb_zip_relrisk(sum(y2), mu2, p2, d2), 2)
})

test_that("score_zip_eb", {
  # This input should give q_hat = 1 immediately
  in1 <- list(y = c(1,0), mu = c(0.5, 2), p = c(0.1, 0.2))
  out1_expected <- list(0, 1) # loglihood score and relative risk estimate
  out1_actual <- score_zip_eb(in1$y, in1$mu, in1$p)
  expect_equal(out1_actual[[1]], out1_expected[[1]])
  expect_equal(out1_actual[[2]], out1_expected[[2]])
  
  # This input holds no zeros and should give q_hat = 2
  in2 <- list(y = c(2, 2), mu = c(1, 1), p = c(0.1, 0.2))
  out2_expected <- list(2 * (zip_lpmf(2, 2, 0.1) - zip_lpmf(2, 1, 0.1)), 2)
  out2_actual <- score_zip_eb(in2$y, in2$mu, in2$p, 1e-3)
  expect_equal(as.numeric(out2_actual[[1]]), as.numeric(out2_expected[[1]]))
  expect_equal(out2_actual[[2]], out2_expected[[2]])
  
})

test_that("calc_all_zip_eb", {
  # Single timepoint
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 2), nrow = 1),
    probs = matrix(c(0.1, 0.2), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))
  
  actual1 <- calc_all_zip_eb(in1$counts, in1$baselines, in1$probs, 
                             in1$zones_flat - 1, in1$zone_lengths)
  expected1_score <- c(zip_lpmf(1, 1, 0.1) - zip_lpmf(1, 0.5, 0.1), 0, 0)
  expect_equal(actual1$score, expected1_score)
  expect_equal(actual1$relrisk, c(2, 1, 1))
  
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
  
  actual2 <- calc_all_zip_eb(in2$counts, in2$baselines, in2$probs, 
                             in2$zones_flat - 1, in2$zone_lengths)
  
  expected2_relrisk <- c(
    2, 
    1, 
    1, 
    3, 
    1,  
    4 / (3 + 2 * (1 - est_zip_zero_indicator(2, 0.2, 1))),
    3 / (1 + 0.5 * (1 - est_zip_zero_indicator(0.5, 0.1, 1))),
    21 / (4 + 2 * (1 - est_zip_zero_indicator(
      2, 0.2, 21 / (4 + 2 * (1 - est_zip_zero_indicator(2, 0.2, 1)))))))
  
  expected2_score <- c(
    # Duration = 1
    zip_lpmf(1, expected2_relrisk[1] * 0.5, 0.1) - 
      zip_lpmf(1, 0.5, 0.1),
     0, 0,
    # Duration = 2
     zip_lpmf(3, expected2_relrisk[4] * 1, 0.1) - 
       zip_lpmf(3, 1, 0.1),
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
       zip_loglihood(c(0, 1, 20), rep(2, 3), rep(0.2, 3), 1))
  
  expect_equal(actual2$score[-9], expected2_score)
  expect_equal(actual2$relrisk[-9], expected2_relrisk)
})

test_that("calc_one_zip_eb", {
  # 1 timepoint
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 2), nrow = 1),
    probs = matrix(c(0.1, 0.2), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))
  
  actual1 <- calc_one_zip_eb(in1$counts, in1$baselines, in1$probs, 
                             in1$zones_flat - 1, in1$zone_lengths)
  expected1_score <- zip_lpmf(1, 1, 0.1) - 
                       zip_lpmf(1, 0.5, 0.1)
  expect_equal(actual1$score, expected1_score)
  expect_equal(actual1$relrisk, 2)
  expect_equal(actual1$zone, 1)
  
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
  
  actual2 <- calc_one_zip_eb(in2$counts, in2$baselines, in2$probs, 
                             in2$zones_flat - 1, in2$zone_lengths)
  
  expected2_relrisk <- 21 / (4 + 2 * (1 - est_zip_zero_indicator(
    2, 0.2, 21 / (4 + 2 * (1 - est_zip_zero_indicator(2, 0.2, 1))))))
  expected2_score <- 
    zip_loglihood(c(0, 1, 20), rep(2, 3), rep(0.2, 3), expected2_relrisk) - 
    zip_loglihood(c(0, 1, 20), rep(2, 3), rep(0.2, 3), 1)
  expect_equal(actual2$score, expected2_score)
  expect_equal(actual2$relrisk, expected2_relrisk)
  expect_equal(actual2$zone, 2)
})