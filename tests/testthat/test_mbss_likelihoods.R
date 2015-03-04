context("MBSS likelihood functions")



# Spatial log-likelihood with uniform prior ---------------------------

test_that("spatial_llh: value of llh in table has changed correctly", {
  lr_input <- data.table(region = rep(c(1,2,1,2), 2),
                         location = rep(c(1,2,3,3), 2),
                         event = rep(1:2, each = 4),
                         llh = 1:8)
  
  # add 1, subtract 2, add 3 = add 2
  spatial_llh(lr_input, 
                      null_llh = 1, 
                      L = exp(2),
                      max_duration = exp(-3))
  
  expect_equal(lr_input[, llh], 
               1:8 + 2)
})


# Space-time log-likelihood with uniform prior ---------------------------

test_that("spacetime_llh: value of llh in table has changed correctly", {
  fullr <- data.table(event = rep(1, 12),
                      region = rep(1:3, each = 4),
                      severity = rep(1:2, 3, each = 2),
                      time = rep(0:1, 6),
                      key = c("event", "region", "time"))
  fullr[, llr := log(rep(1:6, each = 2) / 2)]
  setkeyv(fullr, c("event", "region", "severity"))
  
    
  # add 2, subtract 1 from logsumexp(llh)
  stlh <- spacetime_llh(fullr, 
                        null_llh = 2, 
                        L = exp(1))
  
  expect_equal(stlh[, st_llh], 
               log(1:6) + 2 - 1)
})


# Log-likelihood under null hypothesis of no events ---------------------------

test_that("sums correctly", {
  densities <- data.table(stream = rep(c(1,2,1,2), 2),
                          location = rep(c(1,2,3,3), 2),
                          event = rep(1:2, each = 4),
                          density = 1:8)
  
  expect_equal(null_llh(densities), 
               sum(1:8))
})
