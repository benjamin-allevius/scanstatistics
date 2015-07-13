context("ZIP statistic tests")

test_that("estimate_d", {
  p <- c(2, 3, 4) / 10
  mu <- c(3.5, 5.5, 4.5)
  y <- c(2, 6, 0)
  expected <- c(0, 0, p[3] / (p[3] + (1-p[3]) * exp(-mu[3])))
  actual <- estimate_d(p, mu, y)
  expect_equal(actual, expected)
})

test_that("zip_statistic_term", {
  q <- 1.5
  p <- c(2, 3, 4) / 10
  mu <- c(3.5, 5.5, 4.5)
  y <- c(2, 6, 0)
  dstar <- c(0, 0, 0.99)
  ddagger <- c(0, 0, 0.98)
  actual <- zip_statistic_term(q, p, dstar, ddagger, mu, y) 
  expect_equal(length(actual), 3)
})

test_that("zip_em_estimates", {
  p <- 1:4 / 10
  mu <- c(2, 3.5, 5.5, 4.5)
  y <- c(0, 2, 6, 0)
  actual <- zip_em_estimates(p, mu, y) 
  expect_equal(length(actual), 3)
})

test_that("zip_statistic", {
  p <- 1:4 / 10
  mu <- c(2, 3.5, 5.5, 4.5)
  y <- c(0, 2, 6, 0)
  actual <- zip_em_estimates(p, mu, y) 
  expect_equal(length(actual), 3)
})


test_that("calc_zipstat_over_duration", {
  table <- table_creator(list(location = 1:2, duration = 1:3), 
                         keys = c("location", "duration"))
  table[, mean := 1:6 + 0.5]
  table[, count := c(1,3,2, 7, 3, 10)]
  regions <- sets::set(sets::as.set(1L), 
                       sets::as.set(2L),
                       sets::as.set(1:2))
  a <- region_joiner(table, regions = regions, keys = c("region", "duration"))
  # a[, calc_zipstat_over_duration(.SD, 3), by = .(region)]
  
  p <- 1:4 / 10
  mu <- c(2, 3.5, 5.5, 4.5)
  y <- c(0, 2, 6, 0)
  actual <- zip_em(p, mu, y) 
  expect_equal(length(actual), 3)
})