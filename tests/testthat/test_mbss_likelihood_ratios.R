context("MBSS likelihood ratio functions")


test_that("spatial_llr: sums correctly", {
  y <- data.table(region = rep(1:3, each = 6),
                  event = rep(1:2, 3, each = 3),
                  duration = rep(1:3, 6),
                  llr = log(rep(1:6, 3)))
  setkeyv(y, c("region", "event", "duration"))
  
  
  degl <- matrix(log(1:6), ncol = 2)
  
  sllr <- spatial_llr(y, degl)
  expect_equal(sllr[, llr], rep(log(c(1+4+9, 16+25+36)), 3))
})

test_that("spacetime_llr: sums correctly", {
  y <- data.table(location = c(rep(1:2, each = 6), rep(1:2, 6)),
                  region = c(rep(1:2, each = 6), rep(3, 12)),
                  event = rep(1L, 24),
                  duration = c(rep(1:3, 2, each = 2), rep(1:3, each = 4)),
                  stream = c(rep(1:2, 6), rep(1:2, 3, each = 2)),
                  llr = log(1:24))
  setkeyv(y, c("region", "event", "duration", "stream", "location"))
  
  stllr <- spacetime_llr(y)
  expect_equal(stllr[, exp(llr)], c(3, 10, 21, 15, 34, 57, 58, 132, 222))
})

test_that("event_llr: calculates correctly", {
  llrs <- data.table(location = rep(1:2, each = 8),
                     event = rep(1:2, 4, each = 2),
                     duration = rep(1:2, 2, each = 4),
                     stream = rep(1L, 16),
                     effect_logprob = 16:1,
                     llr = 1:16)
  setkeyv(llrs, c("location", "duration", "event"))
  evllr <- event_llr(llrs)
  expect_equal(evllr[, llr], c(rep(17 + log(2), 8)))
})