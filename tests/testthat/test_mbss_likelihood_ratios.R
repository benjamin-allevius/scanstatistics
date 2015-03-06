context("MBSS likelihood ratio functions")

test_that("full_llr: cumsum over time correct when no missing values", {
  llrs <- data.table(region = rep(1:3, each = 3),
                     event = rep(1, 9),
                     time = rep(0:2, 3),
                     severity = rep(1, 9))
  setkeyv(llrs, c("region", "event", "time", "severity"))
  llrs[, llr := 1:9]
  
  fullr <- full_llr(llrs)
  expect_equal(fullr[, llr],
               c(cumsum(1:3), cumsum(4:6), cumsum(7:9)))
})


test_that("sum_over_streams: sums correctly when no missing values", {
  llrs <- data.table(region = rep(1:3, each = 4),
                     event = rep(1, 12),
                     time = rep(0:1, 3, each = 2),
                     severity = rep(1, 12),
                     stream = rep(1:2, 6))
  setkeyv(llrs, c("region", "event", "time", "severity", "stream"))
  llrs[, llr := rep(1:6, each = 2) / 2]
  
  
  stream_sums <- sum_over_streams(llrs)
  expect_equal(stream_sums[, llr], 1:6)
})


test_that("sum_locations_in_region: sums correctly when no missing values", {
  llrs <- data.table(location = c(rep(1:2, each = 4), rep(1:2, 4)),
                     region = c(rep(1:2, each = 4), rep(3, 8)),
                     event = rep(1, 16),
                     time = c(rep(0:1, 2, each = 2), rep(0:1, each = 4)),
                     severity = rep(1, 16),
                     stream = c(rep(1:2, 4), rep(1:2, 2, each = 2)))
  setkeyv(llrs, c("region", "event", "time", "severity", "stream"))
  llrs[, llr := c(1:8, rep(9:12 / 2, each = 2))]
  
  region_sums <- sum_locations_in_region(llrs)
  expect_equal(region_sums[, llr], 1:12)
})



