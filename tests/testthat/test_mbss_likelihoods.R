context("MBSS likelihood functions")



# Spatial log-likelihood with uniform prior ---------------------------

test_that("spatial_llh: value of llh in table has changed correctly", {
  fullr <- data.table(event = rep(1, 12),
                      region = rep(1:3, each = 4),
                      severity = rep(1:2, 3, each = 2),
                      time = rep(0:1, 6),
                      key = c("event", "region", "time"))
  fullr[, llr := log(rep(1:3, each = 4) / 4)]
  
  # add 2, subtract 1, subtract 1 from logsumexp(llh)
  sllh <- spatial_llh(fullr, 
                      null_llh = 2, 
                      L = exp(1),
                      max_duration = exp(3))
  expect_equal(sllh[, llh], 
               log(1:3) + 2 - 1 - 3)
})


# Space-time log-likelihood with uniform prior ---------------------------

test_that("spacetime_llh: calculated correctly", {
  fullr <- data.table(event = rep(1, 12),
                      region = rep(1:3, each = 4),
                      severity = rep(1:2, 3, each = 2),
                      time = rep(0:1, 6),
                      key = c("event", "region", "time"))
  fullr[, llr := log(rep(1:6, each = 2) / 2)]
  setkeyv(fullr, c("event", "region", "severity"))
    
  # add 2, subtract 1 from logsumexp(llh)
  stllh <- spacetime_llh(fullr, 
                         null_llh = 2, 
                         L = exp(1))
  expect_equal(stllh[, llh], 
               log(1:6) + 2 - 1)
})


# Log-likelihood under null hypothesis of no events ---------------------------

test_that("null_llh: sums correctly", {
  densities <- data.table(stream = rep(c(1,2,1,2), 2),
                          location = rep(c(1,2,3,3), 2),
                          event = rep(1:2, each = 4),
                          density = 1:8)
  
  expect_equal(null_llh(densities), 
               sum(1:8))
})
