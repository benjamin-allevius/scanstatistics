context("MBSS likelihood functions")



# Spatial log-likelihood with uniform prior ---------------------------

test_that("output data.table has correct names", {
  lr_input <- data.table(region = rep(c(1,2,1,2), 2),
                         location = rep(c(1,2,3,3), 2),
                         event = rep(1:2, each = 4),
                         llh = 1:8)
  
  spatial_llh_uniform(lr_input, 
                      null_llh = 1, 
                      L = exp(2),
                      max_duration = exp(-3))
  
  expect_equal(names(lr_input), 
               c("region",
                 "location",                  
                 "event",
                 "llh"))
})

test_that("value of llh in table has changed correctly", {
  lr_input <- data.table(region = rep(c(1,2,1,2), 2),
                         location = rep(c(1,2,3,3), 2),
                         event = rep(1:2, each = 4),
                         llh = 1:8)
  
  # add 1, subtract 2, add 3 = add 2
  spatial_llh_uniform(lr_input, 
                      null_llh = 1, 
                      L = exp(2),
                      max_duration = exp(-3))
  
  expect_equal(lr_input[, llh], 
               1:8 + 2)
})


# Space-time log-likelihood with uniform prior ---------------------------

test_that("output data.table has correct names", {
  stlr_input <- data.table(region = rep(c(1,2,1,2), 4),
                           location = rep(c(1,2,3,3), 4),
                           event = rep(1:2, 2, each = 4),
                           duration = rep(1:2, each = 8),
                           llh = 1:16)
  
  spacetime_llh_uniform(stlr_input, 
                null_llh = 1, 
                L = exp(2))
  
  expect_equal(names(stlr_input), 
               c("region",
                 "location",                  
                 "event",
                 "duration",
                 "llh"))
})

test_that("value of llh in table has changed correctly", {
  stlr_input <- data.table(region = rep(c(1,2,1,2), 4),
                           location = rep(c(1,2,3,3), 4),
                           event = rep(1:2, 2, each = 4),
                           duration = rep(1:2, each = 8),
                           llh = 1:16)
  
  # add one, subtract 2 from llh
  spacetime_llh_uniform(stlr_input, 
                null_llh = 1, 
                L = exp(2))
  
  expect_equal(stlr_input[, llh], 
               1:16 - 1)
})


# Log-likelihood under null hypothesis of no events ---------------------------

test_that("sums correctly", {
  llh_stream <- data.table(stream = 1:3, 
                           llh_stream_value= c(1,10,100),
                           key = "stream")
  llh_rest <- data.table(stream = 1:3, 
                          llh_sum_st = c(-1,-10,-100),
                          key = "stream")
  
  expect_equal(null_llh(llh_stream, llh_rest), 
               0)
})


test_that("sums correctly", {
  llh_mnt <- data.table(stream = rep(1:3, each = 4), 
                        location = rep(1:2, 6),
                        time = rep(0:1, 3, each = 2),
                        llh = 1:12,
                        key = "stream")

  expect_equal(null_llh_rest(llh_mnt)[, llh_sum_st], 
               c(sum(1:4), sum(5:8), sum(9:12)))
})