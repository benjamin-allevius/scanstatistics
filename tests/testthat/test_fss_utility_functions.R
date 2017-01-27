context("Fast Subset Scan: Utility Functions")

test_that("aggregate_per_location: sums properly + dimensions correct", {
  A <- matrix(1:6, 3, 2)
  B <- matrix(-(1:6), 3, 2)
  input <- array(c(A, B), dim = c(3, 2, 2))
  expected <- apply(A, 2, cumsum) + apply(B, 2, cumsum)
  actual <- aggregate_per_location(input)
  expect_identical(dim(actual), dim(expected))
  expect_equal(actual, expected)
})

test_that("aggregate_per_stream: sums properly + dimensions correct", {
  A <- matrix(1:6, 3, 2)
  B <- matrix(-(1:6), 3, 2)
  input <- array(c(A, B), dim = c(3, 2, 2))
  expected <- apply(cbind(apply(A, 1, sum), apply(B, 1, sum)), 2, cumsum)
  actual <- aggregate_per_stream(input)
  expect_identical(dim(actual), dim(expected))
  expect_equal(actual, expected)
})

test_that("apply_rowwise: works + dimensions correct", {
  A <- matrix(1:6, 2, 3)
  expected1 <- matrix(c(1, 4, 9, 
                        2, 6, 12), 
                      2, 3, byrow = TRUE)
  expected2 <- c(9, 12)
  actual1 <- apply_rowwise(A, cumsum)
  actual2 <- apply_rowwise(A, sum)
  expect_equal(actual1, expected1)
  expect_equal(actual2, expected2)
  expect_true(is.matrix(actual1))
  expect_true(is.vector(actual2))
})

test_that("prioritize_cols: no ties", {
  A <- matrix(c(-1, 0, 1,
                0, -1, 1,
                1, 0, -1), 
              3, 3, byrow = TRUE)
  expected <- matrix(c(3, 2, 1,
                       3, 1, 2,
                       1, 2, 3),
                     3, 3, byrow = TRUE)
  actual <- prioritize_cols(A)
  expect_equal(actual, expected)
})

test_that("prioritize_cols: ties", {
  A <- matrix(c(1, 1, 0,
                1, 0, 1,
                0, 1, 1), 
              3, 3, byrow = TRUE)
  expected <- matrix(c(1, 2, 3,
                       1, 3, 2,
                       2, 3, 1),
                     3, 3, byrow = TRUE)
  actual <- prioritize_cols(A)
  expect_equal(actual, expected)
})

test_that("reorder_rows: works", {
  A <- matrix(1:9, 3, 3, byrow = TRUE)
  prios <- matrix(c(2, 3, 1,
                    1, 3, 2,
                    3, 2, 1), 
                  3, 3, byrow = TRUE)
  expected <- matrix(c(2, 3, 1,
                       4, 6, 5,
                       9, 8, 7),
                     3, 3, byrow = TRUE)
  actual <- reorder_rows(A, prios)
  expect_equal(actual, expected)
})

test_that("prioritize_and_execute", {
  A <- matrix(1:9, 3, 3, byrow = TRUE)
  prios <- matrix(c(2, 3, 1,
                    1, 3, 2,
                    3, 2, 1), 
                  3, 3, byrow = TRUE)
  f <- function(x, s = 2) (x + 1) * s
  B <- matrix(c(2, 3, 1,
                4, 6, 5,
                9, 8, 7),
              3, 3, byrow = TRUE)
  expected1 <- (B + 1) * 2
  expected2 <- (B + 1) * 5
  actual1 <- prioritize_and_execute(f, A, prios)
  actual2 <- prioritize_and_execute(f, A, prios, s = 5)
  expect_equal(actual1, expected1)
  expect_equal(actual2, expected2)
})

test_that("sum_over_subset: streams", {
  A <- matrix(1:6, 3, 2)
  B <- matrix(-(1:6), 3, 2)
  arr1 <- array(c(A, B), dim = c(3, 2, 2))
  arr2 <- 2 * arr1
  lst <- list(x = arr1, y = arr2)
  expected1 <- list(x = A + B, y = 2*(A + B))
  expected2 <- list(x = A, y = 2 * A)
  actual1 <- sum_over_subset(lst, 1:2, 3)
  actual2 <- sum_over_subset(lst, 1, 3)
  expect_equal(actual1, expected1)
  expect_equal(actual2, expected2)
})

test_that("sum_over_subset: locations", {
  A <- matrix(1:6, 3, 2)
  B <- matrix(-(1:6), 3, 2)
  arr1 <- array(c(A, B), dim = c(3, 2, 2))
  arr2 <- 2 * arr1
  lst <- list(x = arr1, y = arr2)
  expected1 <- list(x = cbind(A[, 1] + A[, 2], B[, 1] + B[, 2]), 
                    y = 2*cbind(A[, 1] + A[, 2], B[, 1] + B[, 2]))
  expected2 <- list(x = cbind(A[, 1], B[, 1]), 
                    y = 2 * cbind(A[, 1], B[, 1]))
  actual1 <- sum_over_subset(lst, 1:2, 2)
  actual2 <- sum_over_subset(lst, 1, 2)
  expect_equal(actual1, expected1)
  expect_equal(actual2, expected2)
})
