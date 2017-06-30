# context("C++ probability functions")
# 
# 
# # ZIP distribution -------------------------------------------------------------
# 
# test_that("zip_lpmf", {
#   expect_equal(zip_lpmf(0, 6, 0.2),
#                log(0.2 + 0.8 * exp(-6)))
#   expect_equal(zip_loglihood(1, 2, 0.2, 3),
#                log(0.8) + 1 * log(6) - lgamma(1 + 1) - 6)
#   
# })
# 
# test_that("zip_loglihood", {
#   expect_equal(zip_loglihood(c(0, 1), c(3, 3), c(0.2, 0.2), 2),
#                log(0.2 + 0.8 * exp(-2 * 3)) + 
#                  log(0.8) + 1 * log(2 * 3) - lgamma(1 + 1) - 2 * 3)
#   
# })
# 
