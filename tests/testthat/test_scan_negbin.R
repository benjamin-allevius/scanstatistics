context("Negative Binomial Scanstatistics")

### General functions ----------------------------------------------------------

test_that("negbin_mcsim", {
  table <- create_table(list(location = 1:2, duration = 1:3), 
                        keys = c("location", "duration"))
  table[, mu := 1:6 + 0.5]
  table[, theta := c(0.2, 0.5, 2.5, 5, 10, 100)]
  table[, count := c(5, 2, 5, 15, 12, 3)]
  zones <- sets::set(sets::as.set(1L), 
                       sets::as.set(2L),
                       sets::as.set(1:2))
  nsims <- 10
  
  actual_ordinary <- negbin_mcsim(table, zones, nsims, version = "ordinary")
  expect_true(length(actual_ordinary) == nsims)
  expect_true(!any(is.na(actual_ordinary)))
  
  actual_increasing <- negbin_mcsim(table, zones, nsims, version = "increasing")
  expect_true(length(actual_increasing) == nsims)
  expect_true(!any(is.na(actual_increasing)))
})

test_that("negbin_score_terms: calculates correctly", {
  table <- create_table(list(location = 1:2, duration = 1:2), 
                        keys = c("location", "duration"))
  x <- c(1, 3, 5, 7)
  m <- c(0, 3, 1, 10) + 0.5
  s <- 1:4
  table[, count := x]
  table[, mu := m]
  table[, overdispersion := s]
  expected_num <- (x - m) / s
  expected_den <- m / s
  actual <- negbin_score_terms(table)
  expect_equal(actual[, num], expected_num)
  expect_equal(actual[, denom], expected_den)
})

test_that("poisson_score_terms: calculates correctly", {
  table <- create_table(list(location = 1:2, duration = 1:2), 
                        keys = c("location", "duration"))
  x <- c(1, 3, 5, 7)
  m <- c(0, 3, 1, 10) + 0.5
  table[, count := x]
  table[, mu := m]
  expected_num <- x - m
  expected_den <- m
  actual <- poisson_score_terms(table)
  expect_equal(actual[, num], expected_num)
  expect_equal(actual[, denom], expected_den)
})

test_that("score_zone_sums: calculates correctly", {
  table <- create_table(list(location = 1:2, duration = 1:2), 
                        keys = c("location", "duration"))
  table[, num := 1:4]
  table[, denom := 5:8]
  zones <- sets::set(sets::as.set(1L), 
                       sets::as.set(2L),
                       sets::as.set(1:2))
  expected_num <- c(1:4, 1+3, 2+4)
  expected_denom <- c(5:8, 5+7, 6+8)
  actual <- score_zone_sums(table, zones)
  expect_equal(actual[, num], expected_num)
  expect_equal(actual[, denom], expected_denom)
})

### Functions for ordinary outbreak/event/anomaly model ------------------------

test_that("negbin_score: calculated correctly", {
  d <- data.table(zone = rep(1:3, each = 3),
                  duration = rep(1:3, 3))
  d[, num := 1:9 - 5]
  d[, denom := 1:9]
  score_num <- c(cumsum((-4):(-2)),
                 cumsum((-1):1),
                 cumsum(2:4))
  score_denom <- sqrt(c(cumsum(1:3),
                        cumsum(4:6),
                        cumsum(7:9)))
  expected <- score_num / score_denom
  actual <- negbin_score(d)
  expect_equal(actual[, statistic], expected)
})

### Functions for outbreak model -----------------------------------------------

test_that("convolute_numerator: calculated correctly", {
  x <- 5:10
  W <- 1:3
  actual <- convolute_numerator(x, W)
  expected <- c(1*5, 
                2*5 + 1*6, 
                3*5 + 2*6 + 1*7)
  expect_equal(actual, expected)
})

test_that("convolute_denominator: calculated correctly", {
  x <- 5:10
  W <- 1:3
  actual <- convolute_denominator(x, W)
  expected <- c(1^2*5, 
                2^2*5 + 1^2*6, 
                3^2*5 + 2^2*6 + 1^2*7)
  expect_equal(actual, expected)
})

test_that("negbin_increasing_score: calculated correctly", {
  d <- data.table(zone = rep(1:3, each = 3),
                  duration = rep(1:3, 3))
  d[, num := 1:9 - 5]
  d[, denom := 1:9]
  score_num <- c(1*(-4), 2*(-4) + 1*(-3), 3*(-4) + 2*(-3) + 1*(-2),
                 1*(-1), 2*(-1) + 1*0, 3*(-1) + 2*0 + 1*1,
                 1*2, 2*2 + 1*3, 3*2 + 2*3 + 1*4)
  score_denom <- sqrt(c(1, 2^2*1 + 1^2*2, 3^2*1 + 2^2*2 + 1^2*3,
                        4, 2^2*4 + 1^2*5, 3^2*4 + 2^2*5 + 1^2*6,
                        7, 2^2*7 + 1^2*8, 3^2*7 + 2^2*8 + 1^2*9))
  expected <- score_num / score_denom
  actual <- negbin_increasing_score(d)
  expect_equal(actual[, statistic], expected)
})