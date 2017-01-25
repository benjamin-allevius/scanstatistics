context("Fast Subset Scan")

test_that("aggregate_per_location: sums properly + dimensions correct", {
  A <- matrix(1:6, 3, 2)
  B <- matrix(-(1:6), 3, 2)
  input <- array(c(A, B), dim = c(3, 2, 2))
  expected <- apply(A, 2, cumsum) + apply(B, 2, cumsum)
  actual <- aggregate_per_location(input)
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

test_that("prioritize_locations: no ties", {
  A <- matrix(c(-1, 0, 1,
                0, -1, 1,
                1, 0, -1), 
              3, 3, byrow = TRUE)
  expected <- matrix(c(3, 2, 1,
                       3, 1, 2,
                       1, 2, 3),
                     3, 3, byrow = TRUE)
  actual <- prioritize_locations(A)
  expect_equal(actual, expected)
})

test_that("prioritize_locations: ties", {
  A <- matrix(c(1, 1, 0,
                1, 0, 1,
                0, 1, 1), 
              3, 3, byrow = TRUE)
  expected <- matrix(c(1, 2, 3,
                       1, 3, 2,
                       2, 3, 1),
                     3, 3, byrow = TRUE)
  actual <- prioritize_locations(A)
  expect_equal(actual, expected)
})