context("MBSS marginal probability of data")

test_that("logsumexp_llh_over_regions: llh calculated correctly", {
  sllh <- data.table(region = rep(1:3, 2),
                     event = rep(1:2, each = 3),
                     llh = 1:6)
  setkeyv(sllh, c("region", "event"))
  event_llh <- logsumexp_llh_over_regions(sllh)
  lse <- function(x) log(sum(exp(x)))
  expect_equal(event_llh[, llh], 
               c(lse(1:3),
                 lse(4:6)))
})