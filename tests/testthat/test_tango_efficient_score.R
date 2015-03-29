context("Tango 2011 efficient score functions")

test_that("score_numerator: calculated correctly", {
  x <- 5:10
  W <- 1:3
  actual <- score_numerator(x, W)
  expected <- c(1*5, 
                2*5 + 1*6, 
                3*5 + 2*6 + 1*7)
  expect_equal(actual, expected)
})

test_that("score_numerator: calculated correctly", {
  x <- 5:10
  W <- 1:3
  actual <- score_denominator(x, W)
  expected <- c(1^2*5, 
                2^2*5 + 1^2*6, 
                3^2*5 + 2^2*6 + 1^2*7)
  expect_equal(actual, expected)
})

test_that("efficient_score_terms_nbin: calculated correctly", {
  d <- data.table(region = c(rep(1:3, each = 3), rep(3, 3)),
                  location = c(rep(1:2, each = 3), rep(1:2, 3)),
                  duration = c(rep(1:3, 2), rep(1:3, each = 2)))
  d[, count := 1:12]
  d[, baseline := 5]
  d[, overdispersion := 2]
  expected_num <- c((1:3 - 5)/2, 
                    (4:6 - 5)/2, 
                    sum((7:8 - 5)/2),
                    sum((9:10 - 5)/2),
                    sum((11:12 - 5)/2))
  expected_denom <- c(rep(5/2, 6),
                      rep(5, 3))
  actual <- efficient_score_terms_nbin(d)
  expect_equal(actual[, num], expected_num)
  expect_equal(actual[, denom], expected_denom)
})

test_that("efficient_score_terms_binom: calculated correctly", {
  d <- data.table(region = c(rep(1:3, each = 3), rep(3, 3)),
                  location = c(rep(1:2, each = 3), rep(1:2, 3)),
                  duration = c(rep(1:3, 2), rep(1:3, each = 2)))
  d[, count := 1:12]
  d[, baseline := 5]
  expected_num <- c((1:3 - 5), 
                    (4:6 - 5), 
                    sum(7:8 - 5),
                    sum(9:10 - 5),
                    sum(11:12 - 5))
  expected_denom <- c(rep(5, 6),
                      rep(10, 3))
  actual <- efficient_score_terms_binom(d)
  expect_equal(actual[, num], expected_num)
  expect_equal(actual[, denom], expected_denom)
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
  expect_equal(actual[, score], expected)
})

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
  expect_equal(actual[, score], expected)
})