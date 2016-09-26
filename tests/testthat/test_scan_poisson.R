context("Poisson Scan Statistic")


test_that("poisson_mcsim", {
  table <- create_table(list(location = 1:2, duration = 1:3), 
                        keys = c("location", "duration"))
  table[, mu := 1:6 + 0.5]
  table[, count := c(1,3,2, 7, 3, 10)]
  zones <- sets::set(sets::as.set(1L), 
                       sets::as.set(2L),
                       sets::as.set(1:2))
  
  nsims <- 10
  actual <- poisson_mcsim(table, zones, nsims)
  expect_true(!any(actual < 0))
  expect_true(length(actual) == nsims)
})

test_that("poisson_statistic: calculates correctly", {
  table <- create_table(list(zone = 1:3, duration = 1:2), 
                        keys = c("zone", "duration"))
  table[, mu := 1:6 + 0.5]
  table[, count := c(1,3,2, 7, 3, 10)]
  table[, relrisk := max(1, count / mu), by = .(zone, duration)]
  
  actual <- poisson_statistic(table)[, statistic]
  expected <- c(0, 
                log(1.2) * 3 - 0.2 * 2.5,
                0,
                log(7/4.5) * 7 - (7/4.5 - 1) * 4.5,
                0,
                log(10/6.5) * 10 - (10/6.5 - 1) * 6.5)
  expect_equal(actual, expected)
})

test_that("poisson_relrisk: calculates correctly", {
  table <- create_table(list(zone = 1:3, duration = 1:2), 
                        keys = c("zone", "duration"))
  table[, mu := 1:6 + 0.5]
  table[, count := c(1,3,2, 7, 3, 10)]
  
  actual <- poisson_relrisk(table)[, relrisk]
  expected <- c(1, 3/2.5, 1, 7/4.5, 1, 10/6.5)
  expect_equal(actual, expected)
})
