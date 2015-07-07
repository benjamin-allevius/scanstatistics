context("Aggregation functions")

test_that("region_sum: calculates correctly", {
  table <- table_creator(list(location = 1:2, duration = 1:2), 
                         keys = c("location", "duration"))
  table[, count := 1:4]
  table[, mean := 4:1 + 0.5]
  regions <- sets::set(sets::as.set(1L), 
                       sets::as.set(2L),
                       sets::as.set(1:2))
  table <- region_joiner(table, regions, c("region", "duration"))
  
  actual <- region_sum(table, c("count", "mean"))
  expected_count <- c(1, 2, 3, 4, 1+3, 2+4)
  expected_mean <- c(4.5, 3.5, 2.5, 1.5, 4.5+2.5, 3.5+1.5)
  expect_equal(actual[, count], expected_count)
  expect_equal(actual[, mean], expected_mean)
})


test_that("cumsum_duration: calculates correctly", {
  table <- table_creator(list(region = 1:2, duration = 1:3), 
                         keys = c("region", "duration"))
  table[, count := 0:5]
  table[, mean := 1:6 + 0.5]
  
  actual <- cumsum_duration(table, c("count", "mean"), "region")
  expected_count <- c(0,1,3, 3, 7, 12)
  expected_mean <- c(1.5, 4, 7.5, 4.5, 10, 16.5)
  expect_equal(actual[, count], expected_count)
  expect_equal(actual[, mean], expected_mean)
})