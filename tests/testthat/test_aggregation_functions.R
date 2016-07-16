context("Aggregation functions")

test_that("zone_sum: calculates correctly", {
  table <- create_table(list(location = 1:2, duration = 1:2), 
                        keys = c("location", "duration"))
  table[, count := 1:4]
  table[, mean := 4:1 + 0.5]
  zones <- sets::set(sets::as.set(1L), 
                       sets::as.set(2L),
                       sets::as.set(1:2))
  table <- join_zones(table, zones, c("zone", "duration"))
  
  actual <- zone_sum(table, c("count", "mean"))
  expected_count <- c(1, 2, 3, 4, 1+3, 2+4)
  expected_mean <- c(4.5, 3.5, 2.5, 1.5, 4.5+2.5, 3.5+1.5)
  expect_equal(actual[, count], expected_count)
  expect_equal(actual[, mean], expected_mean)
})


test_that("cumsum_duration: calculates correctly", {
  table <- create_table(list(zone = 1:2, duration = 1:3), 
                        keys = c("zone", "duration"))
  table[, count := 0:5]
  table[, mean := 1:6 + 0.5]
  
  actual <- cumsum_duration(table, c("count", "mean"), "zone")
  expected_count <- c(0,1,3, 3, 7, 12)
  expected_mean <- c(1.5, 4, 7.5, 4.5, 10, 16.5)
  expect_equal(actual[, count], expected_count)
  expect_equal(actual[, mean], expected_mean)
})