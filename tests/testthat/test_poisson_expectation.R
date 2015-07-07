context("Expectation-Based Poisson Scan Statistic")


test_that("poisson_statistic: calculates correctly", {
  table <- table_creator(list(region = 1:3, duration = 1:2), 
                         keys = c("region", "duration"))
  table[, mean := 1:6 + 0.5]
  table[, count := c(1,3,2, 7, 3, 10)]
  table[, relrisk := max(1, count / mean), by = .(region, duration)]
  
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
  table <- table_creator(list(region = 1:3, duration = 1:2), 
                         keys = c("region", "duration"))
  table[, mean := 1:6 + 0.5]
  table[, count := c(1,3,2, 7, 3, 10)]
  
  actual <- poisson_relrisk(table)[, relrisk]
  expected <- c(1, 3/2.5, 1, 7/4.5, 1, 10/6.5)
  expect_equal(actual, expected)
})
