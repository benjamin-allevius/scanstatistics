context("MBSS main functions")

# Define some input ----------------------------------------------------------
set.seed(123)
n_times <- 5
times <- seq(as.POSIXct("2015-03-27 11:00:00 CET"), 
             length.out = n_times, by = "15 mins")
locations <- 1:2
streams <- 1:2
events <- 1:2
counts_and_logdensities <- table_creator(list(time = times, 
                                              location = locations,
                                              stream = streams,
                                              event = events))
# Time-dependent rate
null_rate <- counts_and_logdensities[event == 1, 
  5*stream + location + 
    4*sin((lubridate::hour(time) + lubridate::minute(time)/4) / (2*pi*24))]

# Fill with null and event likelihoods
counts_and_logdensities[, 
  count := rep(rpois(length(null_rate), lambda = null_rate), 2)]
counts_and_logdensities[, 
  null_logdensity := dpois(count, lambda = rep(null_rate, 2), log = TRUE)]
counts_and_logdensities[, 
  event_logdensity := dpois(count, lambda = null_rate + rep(1:2, each = 20), 
                            log = TRUE)]

null_prior <- 0.9
event_priors <- c(0.03, 0.07)
duration_condpriors <- matrix(c(1:n_times, rev(1:n_times)) / sum(1:n_times), 
                              ncol = length(events))

regions <- sets::set(sets::as.set(1L), sets::as.set(2L), sets::as.set(1:2))

# Tests ------------------------------------------------------------------------

test_that("MBSS: works for valid input", {
  m <- MBSS(counts_and_logdensities, 
            regions, 
            null_prior, 
            event_priors, 
            duration_condpriors)
  totalprob <- m$null_posterior + m$event_posteriors[, sum(event_posterior)]
  expect_equal(totalprob, 1)
})

test_that("MBSS: works for valid input (uniform duration)", {
  m <- MBSS(counts_and_logdensities, 
            regions, 
            null_prior, 
            event_priors, 
            "uniform")
  totalprob <- m$null_posterior + m$event_posteriors[, sum(event_posterior)]
  expect_equal(totalprob, 1)
})

test_that("MBSS: throws error when wrong column name", {
  cal <- copy(counts_and_logdensities)
  setnames(cal, "time", "invalid")
  expect_error(MBSS(cal, 
                    regions, 
                    null_prior, 
                    event_priors, 
                    "uniform"))
  
})

test_that("MBSS: throws error when wrong column name", {
  cal <- copy(counts_and_logdensities)
  setnames(cal, "time", "invalid")
  expect_error(MBSS(cal, 
                    regions, 
                    null_prior, 
                    event_priors, 
                    "uniform"))
  
})

test_that("MBSS: throws error when NA present", {
  cal <- copy(counts_and_logdensities)
  cal[1, ]
  expect_error(MBSS(cal, 
                    regions, 
                    null_prior, 
                    event_priors, 
                    "uniform"))
  
})