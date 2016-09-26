context("ZIP statistic tests")

test_that("estimate_d", {
  p <- c(2, 3, 4) / 10
  mu <- c(3.5, 5.5, 4.5)
  y <- c(2, 6, 0)
  expected <- c(0, 0, p[3] / (p[3] + (1-p[3]) * exp(-mu[3])))
  actual <- estimate_d(p, mu, y)
  expect_equal(actual, expected)
})

test_that("zip_statistic_term: good input", {
  q <- 1.5
  p <- c(2, 3, 4) / 10
  mu <- c(3.5, 5.5, 4.5)
  y <- c(2, 6, 0)
  dstar <- c(0, 0, 0.99)
  ddagger <- c(0, 0, 0.98)
  actual <- zip_statistic_term(q, p, dstar, ddagger, mu, y)
  expect_equal(length(actual), 3)
})

test_that("zip_statistic_term: bad input", {
  q <- 1
  p <- c(0, 3, 0) / 10
  mu <- c(1, 5.5, 4.5)
  y <- c(2, 6, 3)
  dstar <- c(0, 0, 0.99)
  ddagger <- c(0, 0, 0.98)
  actual <- zip_statistic_term(q, p, dstar, ddagger, mu, y) 
  expect_equal(length(actual), 3)
})

test_that("zip_em_estimates: no-outbreak input", {
  p <- 1:4 / 10
  mu <- c(2, 3.5, 5.5, 4.5)
  y <- c(0, 2, 6, 0)
  actual <- zip_em_estimates(p, mu, y) # gamlss.dist::dZIP(y, mu, p)
  expect_true(actual$q == 1)
  expect_true(all(actual$dstar[2:3] == 0))
})

test_that("zip_em_estimates: outbreak input", {
  p <- 1:4 / 10
  mu <- c(2, 3.5, 5.5, 4.5)
  y <- c(0, 10, 16, 0)
  eps <- 0.01
  actual <- zip_em_estimates(p, mu, y) # gamlss.dist::dZIP(y, mu, p)
  expect_true(actual$q > 1 + eps)
  expect_true(all(actual$dstar[2:3] == 0))
})

test_that("zip_em_estimates: outbreak input => structural zeros more likely", {
  # Structural zeros should be more likely since counts are generated from ZIP
  # distribution with higher Poisson mu.
  p <- 1:4 / 10
  mu <- c(2, 3.5, 5.5, 4.5)
  actual_no <- zip_em_estimates(p, mu, c(0, 2, 6, 0))
  actual_yes <- zip_em_estimates(p, mu, c(0, 10, 16, 0))
  expect_true(actual_yes$q > actual_no$q)
  expect_true(all(actual_yes$dstar[c(1,4)] > actual_no$dstar[c(1,4)]))
})

test_that("window_zip_statistic", {
  p <- 1:4 / 10
  mu <- c(2, 3.5, 5.5, 4.5)
  y <- c(0, 2, 6, 0)
  actual <- window_zip_statistic(p, mu, y) 
  expect_equal(length(actual), 1)
})


test_that("calc_zipstat_over_duration: works", {
  table <- create_table(list(location = 1:2, duration = 1:3), 
                        keys = c("location", "duration"))
  table[, mu := 1:6 + 0.5]
  table[, p := 1:6 / 20]
  # Counts should correspond to outbreak with duration 2 at location 1
  table[, count := c(5, 10, 1, 4, 5, 0)]
  actual <- calc_zipstat_over_duration(table, 3)
  expected <- c(
    table[duration <= 1L, window_zip_statistic(p, mu, count)],
    table[duration <= 2L, window_zip_statistic(p, mu, count)],
    table[duration <= 3L, window_zip_statistic(p, mu, count)]
  )
  expect_equal(actual$statistic, expected)
})


test_that("zip_statistic: works", {
  table <- create_table(list(location = 1:2, duration = 1:3), 
                        keys = c("location", "duration"))
  table[, mu := 1:6 + 0.5]
  table[, p := 1:6 / 20]
  # Counts should correspond to outbreak with duration 2 at location 1
  table[, count := c(5, 10, 1, 4, 5, 0)] 
  # table[, gamlss.dist::dZIP(count, mu, p)]
  zones <- sets::set(sets::as.set(1L), 
                       sets::as.set(2L),
                       sets::as.set(1:2))
  dt <- join_zones(table, zones = zones, keys = c("zone", "duration"))
  actual <- zip_statistic(dt, maxdur = 3)
  # actual <- dt[, calc_zipstat_over_duration(.SD, 3), by = .(zone)]
  expected <- c(
    dt[1L, window_zip_statistic(p, mu, count)],
    dt[1:2, window_zip_statistic(p, mu, count)],
    dt[1:3, window_zip_statistic(p, mu, count)],
    dt[4L, window_zip_statistic(p, mu, count)],
    dt[4:5, window_zip_statistic(p, mu, count)],
    dt[4:6, window_zip_statistic(p, mu, count)],
    dt[7:8, window_zip_statistic(p, mu, count)],
    dt[7:10, window_zip_statistic(p, mu, count)],
    dt[7:12, window_zip_statistic(p, mu, count)]
  )
  expect_equal(actual[, statistic], expected)
})

test_that("zip_mcsim", {
  table <- create_table(list(location = 1:2, duration = 1:3), 
                        keys = c("location", "duration"))
  table[, mu := 1:6 + 0.5]
  table[, p := 1:6 / 20]
  zones <- sets::set(sets::as.set(1L), 
                       sets::as.set(2L),
                       sets::as.set(1:2))
  nsims <- 10
  # set.seed(25)
  actual <- zip_mcsim(table, zones, nsims, maxdur = 3)
  expect_true(length(actual) == nsims)
  # expect_true(!any(actual < 0))
})

# test_that("zip_scanstatistic", {
#   table <- create_table(list(location = 1:2, duration = 1:3), 
#                         keys = c("location", "duration"))
#   table[, mu := 1:6 + 0.5]
#   table[, p := 1:6 / 20]
#   # Counts should correspond to outbreak with duration 2 at location 1
#   table[, count := c(5, 10, 1, 4, 5, 0)] 
#   # table[, gamlss.dist::dZIP(count, mu, p)]
#   zones <- sets::set(sets::as.set(1L), 
#                        sets::as.set(2L),
#                        sets::as.set(1:2))
#   nsims <- 9
#   actual <- zip_scanstatistic(table, zones, nsims)
#   expect_equal(actual$mlc[, zone], 1L)
#   expect_equal(actual$mlc[, duration], 2L)
#   })