context("Fast Subset Scan: Exponential Family Functions")

test_that("poisson_qmax", {
  expect_equal(poisson_qmax(5, 5), 1)
  # From Wolfram Alpha: solve 10*log(x)+5*(1-x)=0
  expect_equal(round(poisson_qmax(10, 5), 5), 3.51286)
})

test_that("poisson_priority", {
  B <- matrix(c(1, 3,
                2, 4), 2, 2, byrow = TRUE)
  C <- matrix(c(2, 3,
                1, 5), 2, 2, byrow = TRUE)
  expected <- matrix(c(poisson_qmax(2, 1), 1,
                       1,                  poisson_qmax(8, 7)),
                     2, 2, byrow = TRUE)
  actual <- poisson_priority(C, B, poisson_qmax)
  expect_equal(actual, expected)
})

test_that("poisson_score", {
  pri_mat <- matrix(c(1, 2,
                      2, 1), 
                    2, 2, byrow = TRUE)
  B <- matrix(c(1, 3,
                2, 4), 2, 2, byrow = TRUE)
  C <- matrix(c(2, 3,
                1, 5), 2, 2, byrow = TRUE)
  expected <- matrix(c(poisson_lambda(2, 1), poisson_lambda(2+3,1+3),
                       poisson_lambda(3+5, 3+4), poisson_lambda(2+1+3+5,
                                                                1+2+3+4)),
                     2, 2, byrow = TRUE)
  actual <- poisson_score(C, B, pri_mat)
  expect_equal(actual, expected)
})


