context("utility functions")

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