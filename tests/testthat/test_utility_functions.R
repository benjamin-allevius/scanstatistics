context("Utility functions")

# Functions to do with data.table keys -----------------------------------------
test_that("key is returned if present", {
  DT <- data.table(x = 1:3, y = 1:3, key = "x")
  expect_equal(getkeys(DT), "x")
})

test_that("keys are returned if present", {
  DT <- data.table(x = 1:3, y = 1:3, key = c("x", "y"))
  expect_equal(getkeys(DT), c("x", "y"))
})

test_that("no keys outputs NULL", {
  DT <- data.table(x = 1:3, y = 1:3)
  expect_null(getkeys(DT))
})


test_that("first key correct", {
  DT <- data.table(x = 1:3, y = 1:3, key = c("x", "y"))
  expect_true(first_keys_are_equal(DT, "x"))
})

test_that("first few keys correct", {
  DT <- data.table(x = 1:3, y = 1:3, key = c("x", "y"))
  expect_true(first_keys_are_equal(DT, "x"))
})

test_that("all keys correct", {
  DT <- data.table(x = 1:3, y = 1:3, z = 1:3, key = c("x", "y"))
  expect_true(first_keys_are_equal(DT, c("x", "y")))
})

test_that("FALSE if no keys", {
  DT <- data.table(x = 1:3, y = 1:3)
  expect_false(first_keys_are_equal(DT, c("x", "y")))
})

test_that("FALSE if checking too many keys", {
  DT <- data.table(x = 1:3, y = 1:3, key = c("x", "y"))
  expect_false(first_keys_are_equal(DT, c("x", "y", "z")))
})

# Convenience ------------------------------------------------------------------
test_that("add_llr: works correctly", {
  n_times <- 2
  times <- seq(as.POSIXct("2015-03-27 11:00:00 CET"), 
               length.out = n_times, by = "15 mins")
  locations <- 1:2
  streams <- 1:2
  events <- 0:2
  loglikelihoods <- table_creator(list(time = times, 
                                       location = locations,
                                       stream = streams,
                                       event = events))
  setkeyv(loglikelihoods, c("event", "location", "stream", "time"))
  loglikelihoods[, loglikelihood := 1:24]
  expected <- loglikelihoods[event != 0][, 
    llr := loglikelihood - rep(1:8, 2)]
  setkeyv(expected, c("location", "event", "stream", "time"))
  expect_equal(add_llr(loglikelihoods, 0)[, llr], expected[, llr])
})