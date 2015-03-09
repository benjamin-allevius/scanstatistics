context("MBSS marginal probability of data")

test_that("marginal_data_logprobability: llh calculated correctly", {
  expect_equal(marginal_data_logprobability(log(2), log(3), log(4)),
               log(10))
})

test_that("data_logprob_if_event: llh calculated correctly", {
  event_llh <- data.table(event = 1:3, llh = 1:3 / 2, key = "event")
  
  dataprob_event <- data_logprob_if_event(event_llh, 3:1 / 2, exp(1))
  expect_equal(dataprob_event, 5)
})


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