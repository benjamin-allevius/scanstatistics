context("EB Poisson statistic tests")

test_that("calc_all_poisson_eb", {
  # Single timepoint
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 2), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))
  
  actual1 <- calc_all_poisson_eb(
    in1$counts, in1$baselines, in1$zones_flat - 1, in1$zone_lengths)
  expected1_score <- c(poisson_lpmf(1, 1) - poisson_lpmf(1, 0.5), 0, 0)
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
  
  actual2 <- calc_all_poisson_eb(in2$counts,
                                 apply(in2$baselines, 2, cumsum),
                                 in2$zones_flat - 1, in2$zone_lengths)
  
  expected2_relrisk <- c(1 / 0.5, 1, 1,
                         3, 1, 1,
                         3 / 1.5, 21 / 6, 24 / 7.5)
  
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
  
  expect_equal(actual2$score, expected2_score)
  expect_equal(actual2$relrisk, expected2_relrisk)
})
