context("MBSS likelihood ratio functions")

# Log-likelihood ratios: event of type k and severity l vs no event.
# For all locations and times, log-likelihood ratios are specified in terms 
# of one part that only depends on the data stream m, 
# and event type and severity through impact factor x;
# and another part which also depends on the location i.


test_that("sum_lr_over_time: calculated correctly when no missing values", {
  fullr <- data.table(time = rep(0:3, 3),
                      event = rep(1, 12),
                      severity = rep(1, 12),
                      region = rep(1:3, each = 4), 
                      lq = log(rep(1:4, 3)))
  setkeyv(fullr, c("time", "event", "severity", "region"))
  res <- sum_lr_over_time(fullr)
  expect_equal(res[, lr_timesum],
               rep(1 * (1 + 2 * (1 + 3 * (1 + 4))), 
                   3))
})


test_that("full_llr: cumsum over time correct when no missing values", {
  d1 <- data.table(time = rep(0:3, 3),
                   event = rep(1, 12),
                   severity = rep(1, 12),
                   region = rep(1:3, each = 4),
                   lq = c(0:3, 10:13, 100:103),
                   key = c("time", "event", "severity", "region"))
  
  duration_llr <- full_llr(d1)
  expect_equal(duration_llr[, llr],
               c(cumsum(0:3), cumsum(10:13), cumsum(100:103)))
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



