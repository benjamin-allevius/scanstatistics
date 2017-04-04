context("ZIP statistic tests")

test_that("incomplete_loglihood_term", {
  expect_equal(incomplete_loglihood_term(0, 3, 0.2, 2),
               log(0.2 + 0.8 * exp(-2 * 3)))
  expect_equal(incomplete_loglihood_term(1, 3, 0.2, 2),
               log(0.8) + 1 * log(2 * 3) - lgamma(1 + 1) - 2 * 3)

})

test_that("incomplete_loglihood", {
  expect_equal(incomplete_loglihood(c(0, 1), c(3, 3), c(0.2, 0.2), 2),
               log(0.2 + 0.8 * exp(-2 * 3)) + 
               log(0.8) + 1 * log(2 * 3) - lgamma(1 + 1) - 2 * 3)
  
})

test_that("estimate_q", {
  y1 <- c(1, 0); mu1 <- c(0.5, 2); p1 <- c(0.1, 0.2); d1 <- c(0, 0.65)
  expect_equal(estimate_q(sum(y1), mu1, p1, d1), 1)
  
  y2 <- c(2, 2); mu2 <- c(1, 1); p2 <- c(0.1, 0.2); d2 <- c(0, 0)
  expect_equal(estimate_q(sum(y2), mu2, p2, d2), 2)
})

test_that("zip_em_algo", {
  # This input should give q_hat = 1 immediately
  in1 <- list(y = c(1,0), mu = c(0.5, 2), p = c(0.1, 0.2))
  out1_expected <- list(0, 1) # loglihood score and relative risk estimate
  out1_actual <- zip_em_algo(in1$y, in1$mu, in1$p)
  expect_equal(out1_actual[[1]], out1_expected[[1]])
  expect_equal(out1_actual[[2]], out1_expected[[2]])
  
  # This input holds no zeros and should give q_hat = 2
  in2 <- list(y = c(2,2), mu = c(1, 1), p = c(0.1, 0.2))
  out2_expected <- list(2 * (incomplete_loglihood_term(2, 1, 0.1, 2) - 
                               incomplete_loglihood_term(2, 1, 0.1, 1)), 2) 
  out2_actual <- zip_em_algo(in2$y, in2$mu, in2$p, 1e-3)
  expect_equal(out2_actual[[1]], out2_expected[[1]])
  expect_equal(out2_actual[[2]], out2_expected[[2]])
  
})