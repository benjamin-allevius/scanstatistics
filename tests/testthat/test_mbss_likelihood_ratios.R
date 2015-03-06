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
  kc <- c("stream", "time", "event", "severity")
  
  # Part of LLR dependent on region
  d1 <- data.table(stream = rep(1:2, 6),
                   time = rep(0:1, 3, each = 2),
                   event = rep(1, 12),
                   severity = rep(1, 12),
                   region =  rep(1:3, each = 4))
  setkeyv(d1, kc)
  d1[, slg := rep(1:6, 2) / 2]
  
  # Part of LLR independent of region
  d2 <- data.table(stream = rep(1:2, each = 2),
                   time = rep(0:1, 2),
                   event = rep(1, 4),
                   severity = rep(1, 4),
                   key = kc)
  d2[, lf := c(-2, -1, 2, 1)]
  
  stream_term <- sum_over_streams(d1, d2)
  expect_equal(stream_term[, lq], 1:6)
})


test_that("sum_locations_in_region: sums correctly when no missing values", {
  tab <- data.table(location = rep(1:2, each = 8),
                    region = rep(c(1,3,2,3), each = 4),
                    stream = rep(1:2, 8),
                    time = rep(0:1, 4, each = 2),
                    severity = rep(1, 16),
                    event = rep(1, 16))
  kc <- c("stream", "time", "severity", "event")
  setkeyv(tab, c("region", kc))
  
  tab[, lg := c(1:8, rep(c(10, 100, 1000, 10000) / 2, each = 2))]
  region_sums <- sum_locations_in_region(tab)
  
  # Sort by the sums, for easier comparison
  setkey(region_sums, "slg")
  expect_equal(region_sums[, slg], 
               c(1:8, c(10, 100, 1000, 10000)))
})



