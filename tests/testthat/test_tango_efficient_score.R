context("Tango 2011 efficient score functions")

### General functions ----------------------------------------------------------

test_that("efficient_score_terms_nbin: calculates correctly", {
  table <- table_creator(list(location = 1:2, duration = 1:2), 
                         keys = c("location", "duration"))
  x <- c(1, 3, 5, 7)
  m <- c(0, 3, 1, 10) + 0.5
  s <- 1:4
  table[, count := x]
  table[, mean := m]
  table[, overdispersion := s]
  expected_num <- (x - m) / s
  expected_den <- m / s
  actual <- efficient_score_terms_nbin(table)
  expect_equal(actual[, num], expected_num)
  expect_equal(actual[, denom], expected_den)
})

test_that("efficient_score_terms_poisson: calculates correctly", {
  table <- table_creator(list(location = 1:2, duration = 1:2), 
                         keys = c("location", "duration"))
  x <- c(1, 3, 5, 7)
  m <- c(0, 3, 1, 10) + 0.5
  table[, count := x]
  table[, mean := m]
  expected_num <- x - m
  expected_den <- m
  actual <- efficient_score_terms_poisson(table)
  expect_equal(actual[, num], expected_num)
  expect_equal(actual[, denom], expected_den)
})

test_that("efficient_score_region_sums: calculates correctly", {
  table <- table_creator(list(location = 1:2, duration = 1:2), 
                         keys = c("location", "duration"))
  table[, num := 1:4]
  table[, denom := 5:8]
  regions <- sets::set(sets::as.set(1L), 
                       sets::as.set(2L),
                       sets::as.set(1:2))
  expected_num <- c(1:4, 1+3, 2+4)
  expected_denom <- c(5:8, 5+7, 6+8)
  actual <- efficient_score_region_sums(table, regions)
  expect_equal(actual[, num], expected_num)
  expect_equal(actual[, denom], expected_denom)
})

### Functions for hotspot model ------------------------------------------------

test_that("hotspot_efficient_score: calculated correctly", {
  d <- data.table(region = rep(1:3, each = 3),
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
  actual <- hotspot_efficient_score(d)
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

test_that("outbreak_efficient_score: calculated correctly", {
  d <- data.table(region = rep(1:3, each = 3),
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
  actual <- outbreak_efficient_score(d)
  expect_equal(actual[, statistic], expected)
})