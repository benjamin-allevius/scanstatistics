context("MBSS likelihood functions")


# Space-time log-likelihood with uniform prior ---------------------------

test_that("output data.table has correct names", {
  stlr_input <- data.table(region = rep(c(1,2,1,2), 2),
                           location = rep(c(1,2,3,3), 2),
                           event = rep(1:2, each = 4),
                           llh = c(1:4, 11:14))
  
  spacetime_llh_uniform(stlr_input, 
                null_llh = 1, 
                n_severities = exp(2))
  
  expect_equal(names(stlr_input), 
               c("region",
                 "location",                  
                 "event",
                 "llh"))
})

test_that("table has changed correctly", {
  stlr_input <- data.table(region = rep(c(1,2,1,2), 2),
                           location = rep(c(1,2,3,3), 2),
                           event = rep(1:2, each = 4),
                           llh = c(1:4, 11:14))
  
  # add one, subtract 2 from llh
  spacetime_llh_uniform(stlr_input, 
                null_llh = 1, 
                n_severities = exp(2))
  
  expect_equal(stlr_input[, llh], 
               c(1:4, 11:14) - 1)
})


# Spatial log-likelihood with uniform prior ---------------------------

# Input for space-time log-likelihood function
stlr_input <- data.table(region = rep(c(1,2,1,2), 2),
                         location = rep(c(1,2,3,3), 2),
                         event = rep(1:2, each = 4),
                         llh = c(1:4, 11:14))