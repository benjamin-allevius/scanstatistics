context("MBSS main functions")

# Define some input ----------------------------------------------------------
set.seed(123)
n_times <- 5
times <- seq(as.POSIXct("2015-03-27 11:00:00 CET"), 
             length.out = n_times, by = "15 mins")
locations <- 1:2
streams <- 1:2
events <- 0:2
loglikelihoods <- table_creator(list(time = times, 
                                     location = locations,
                                     stream = streams,
                                     event = events))
# Time-dependent rate
null_rate <- loglikelihoods[event == 0, 
  5*stream + location + 
    8*sin((lubridate::hour(time) + lubridate::minute(time)/4) / (2*pi*24))]

# Fill with null and event likelihoods
counts <- rpois(length(null_rate), lambda = null_rate)
counts[(length(counts) - 2):length(counts)] <- 22

# Fill with null and event likelihoods
loglikelihoods[, loglikelihood := NaN]

# Null log-likelihoods
loglikelihoods[event == 0, 
               loglikelihood := dpois(counts, lambda = null_rate, log = TRUE)]

# Event log-likelihoods
loglikelihoods[event == 1, 
               loglikelihood := dpois(counts, lambda = null_rate + 5, log = TRUE)]
loglikelihoods[event == 2, 
               loglikelihood := dpois(counts, lambda = null_rate + 10, log = TRUE)]

null_prior <- 0.9
event_priors <- c(0.02, 0.05)
duration_condpriors <- matrix(c(1:n_times, rev(1:n_times)) / sum(1:n_times), 
                              ncol = length(events) - 1)

regions <- sets::set(sets::as.set(1L), sets::as.set(2L), sets::as.set(1:2))

# Tests: valid input -----------------------------------------------------------

test_that("MBSS: works for valid input", {
  m <- MBSS(loglikelihoods, 
            regions, 
            null_prior, 
            event_priors, 
            duration_condpriors)
  totalprob <- m$null_posterior + m$event_posteriors[, sum(event_posterior)]
  expect_equal(totalprob, 1)
})

test_that("MBSS: works for valid input (uniform duration)", {
  m <- MBSS(loglikelihoods, 
            regions, 
            null_prior, 
            event_priors, 
            "uniform")
  totalprob <- m$null_posterior + m$event_posteriors[, sum(event_posterior)]
  expect_equal(totalprob, 1)
})

test_that("MBSS: works for valid input; events as strings", {
  cal <- copy(loglikelihoods)
  cal[, event := rep(c("null", "e1", "e2"), each = 20)]
  event_priors <- c(e1 = 2, e2 = 5)
  duration_condpriors <- as.data.frame(duration_condpriors)
  names(duration_condpriors) <- names(event_priors)
  m <- MBSS(cal, 
            regions, 
            null_prior, 
            event_priors, 
            duration_condpriors)
  totalprob <- m$null_posterior + m$event_posteriors[, sum(event_posterior)]
  expect_equal(totalprob, 1)
})

test_that("MBSS: works for valid input; events as strings, uniform duration", {
  cal <- copy(loglikelihoods)
  cal[, event := rep(c("null", "e1", "e2"), each = 20)]
  event_priors <- data.frame(e1 = 2, e2 = 5)
  m <- MBSS(cal, 
            regions, 
            null_prior, 
            event_priors, 
            "uniform")
  totalprob <- m$null_posterior + m$event_posteriors[, sum(event_posterior)]
  expect_equal(totalprob, 1)
})

# Tests: invalid input ---------------------------------------------------------

test_that("MBSS: throws error when input not data.table", {
  cal <- copy(loglikelihoods)
  class(cal) <- "data.frame"
  expect_error(MBSS(cal, regions, null_prior, event_priors, "uniform"),
               regexp = "The log-likelihoods must be supplied as a data.table.")
})

test_that("MBSS: throws error when wrong column name", {
  cal <- copy(loglikelihoods)
  setnames(cal, "time", "invalid")
  expect_error(MBSS(cal, regions, null_prior, event_priors, "uniform"),
               regexp = "The data.table containing")
})

test_that("MBSS: detects positive log-likelihood", {
  cal <- copy(loglikelihoods)
  cal[1, loglikelihood := 2]
  expect_error(MBSS(cal, regions, null_prior, event_priors, "uniform"),
               regexp = "Log-likelihoods should be")
})

test_that("MBSS: throws error when NA present", {
  cal <- copy(loglikelihoods)
  cal[1, null_loglikelihood := NA]
  expect_error(MBSS(cal, regions, null_prior, event_priors, "uniform"),
               regexp = "No support for missing values yet.")
})

test_that("MBSS: throws error when null_prior incorrectly specified", {
  expect_error(MBSS(loglikelihoods, regions, 2, event_priors, "uniform"),
               regexp = "Prior null hypothesis probability")
  expect_error(MBSS(loglikelihoods, regions, -1, event_priors, "uniform"),
               regexp = "Prior null hypothesis probability")
  expect_error(MBSS(loglikelihoods, regions, "2", event_priors, "uniform"),
               regexp = "Prior null hypothesis probability")
  expect_error(MBSS(loglikelihoods, regions, c(1,1)/2, event_priors, "uniform"),
               regexp = "Prior null hypothesis probability")
})