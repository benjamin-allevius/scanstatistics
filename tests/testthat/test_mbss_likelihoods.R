context("MBSS likelihood functions")



# Spatial log-likelihood with uniform prior ---------------------------

test_that("spatial_llh_uniform: llh calculated correctly", {
  fullr <- data.table(event = rep(1, 12),
                      region = rep(1:3, each = 4),
                      severity = rep(1:2, 3, each = 2),
                      duration = rep(1:2, 6),
                      key = c("event", "region", "duration"))
  fullr[, llr := log(rep(1:3, each = 4) / 4)]
  
  # add 2, subtract 1, subtract 1 from logsumexp(llh)
  sllh <- spatial_llh_uniform(fullr, 
                              null_llh = 2, 
                              L = exp(1),
                              max_duration = exp(3))
  expect_equal(sllh[, llh], 
               log(1:3) + 2 - 1 - 3)
})


# Space-duration log-likelihood with uniform prior ---------------------------

test_that("spacetime_llh_uniform: calculated correctly", {
  fullr <- data.table(region = rep(1:3, each = 4),
                      event = rep(1, 12),
                      severity = rep(1:2, 3, each = 2),
                      duration = rep(1:2, 6))
  setkeyv(fullr, c("region", "event", "severity"))
  fullr[, llr := log(c(1,2,1,2,3,4,3,4,5,6,5,6) / 2)]
      
  # add 2, subtract 1 from logsumexp(llh)
  stllh <- spacetime_llh_uniform(fullr, 
                                 null_llh = 2, 
                                 L = exp(1))
  expect_equal(stllh[, llh], 
               log(1:6) + 2 - 1)
})


# Log-likelihood under null hypothesis of no events ---------------------------

test_that("null_llh: calculated correctly", {
  densities <- data.table(stream = rep(c(1,2,1,2), 2),
                          location = rep(c(1,2,3,3), 2),
                          event = rep(1:2, each = 4),
                          density = 1:8)
  
  expect_equal(null_llh(densities), 
               sum(1:8))
})
