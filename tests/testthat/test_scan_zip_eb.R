context("ZIP statistic tests")

test_that("incomplete_loglihood_term", {
  expect_equal(incomplete_loglihood_term(0, 3, 0.2, 2),
               log(0.2 + 0.8 * exp(-2 * 3)))
  expect_equal(incomplete_loglihood_term(1, 3, 0.2, 2),
               log(0.8) + 1 * log(2 * 3) - lgamma(1 + 1) - 2 * 3)

})