context("MBSS posterior probabilities")

test_that("spacetime_logposterior: calculated correctly", {
  st_lps <- data.table(region = rep(1:3, each = 4),
                       event = rep(1:2, 3, each = 2),
                       time = rep(0:1, 6), 
                       llh = rep(c(0.5, 2), 3, each = 2))
  st_lps[, llh := llh + 1:12]
  setkeyv(st_lps, c("region", "event", "time"))
  
  spacetime_logposterior(st_lps, c(0.5, -1), exp(1), exp(2), 3)
  expect_equal(st_lps[, posterior_logprob], 1:12 - 5)
})