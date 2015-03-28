context("MBSS posterior probabilities")


test_that("posterior_event_probabilities: calculated correctly", {
  jd <- data.table(event = rep(1:2, each = 2),
                   duration = rep(1:2, 2),
                   duration_event_posterior = 1:4)
  expected <- c(c(1, 2) / 3, c(3, 4) / 7)
  actual <- posterior_duration_givn_event(jd)[, duration_condposterior]
  expect_equal(actual, expected)
})

test_that("posterior_duration_event_jdist: calculated correctly", {
  st <- data.table(region = rep(1:2, each = 4),
                                event = rep(1:2, 2, each = 2),
                                duration = rep(1:2, 4))
  setkeyv(st, c("event", "duration"))
  st[, posterior_logprob := log(1:8)]
  
  actual <- posterior_duration_event_jdist(st)[, duration_event_posterior]
  expected <- c(1+2, 3+4, 5+6, 7+8)
  expect_equal(actual, expected)
})

test_that("posterior_event_probabilities: calculated correctly", {
  pl <- data.table(event = rep(1:2, each = 2),
                   duration = rep(1:2, 2),
                   duration_event_posterior = 1:4)
  actual <- posterior_event_probabilities(pl)[, event_posterior]
  expected <- c(1+2, 3+4)
  expect_equal(actual, expected)
})

test_that("data_to_nulldata_logratio: calculated correctly", {
  sllr <- data.table(region = rep(1:3, each = 2),
                     event = rep(1:2, 3), 
                     llr = 1:6)
  pH0 <- 2
  n_regions <- 3
  event_logpriors <- c(-1, 1)
  s1 <- log(sum(exp(c(1, 3, 5))))
  s2 <- log(sum(exp(c(2, 4, 6))))
  e1 <- event_logpriors[1] + s1
  e2 <- event_logpriors[2] + s2
  expected <- log(pH0 + exp(-log(n_regions) + log(sum(exp(c(e1, e2))))))
  actual <- data_to_nulldata_logratio(sllr, event_logpriors, pH0, n_regions)
  expect_equal(actual, expected)
})

test_that("spatial_logposterior: calculated correctly", {
  lps <- data.table(region = rep(1:3, each = 2),
                    event = rep(1:2, 3),
                    llr = rep(c(0.5, 2), 3))
  lps[, llr := llr + 1:6]
  setkeyv(lps, c("region", "event"))
  spatial_logposterior(lps, c(0.5, -1), exp(1), 1)
  expect_equal(lps[, posterior_logprob], 1:6 - 1)
})

test_that("spacetime_logposterior: calculated correctly", {
  stllrs <- data.table(region = rep(1:2, each = 4),
                       event = rep(1:2, 2, each = 2),
                       duration = rep(1:2, 4), 
                       llr = 1:8)
  setkeyv(stllrs, c("region", "event", "duration"))
  
  event_logpriors <- c(-1, 1)
  degl <- matrix(1:4, ncol = 2)
  expn_reg <- exp(1)
  data_logratio <- 1
  
  expected <- c(-1, 1, 5, 7, 3, 5, 9, 11)
  actual <- spacetime_logposterior(stllrs, 
                                   event_logpriors, 
                                   degl, 
                                   expn_reg, 
                                   data_logratio)
  expect_equal(actual[, posterior_logprob], expected)
})