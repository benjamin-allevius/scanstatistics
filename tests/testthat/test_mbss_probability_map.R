context("MBSS probability map")

test_that("event_probability_map: calculated correctly", {
  lps <- data.table(region = rep(1:3, each = 2),
                    event = rep(1:2, 3),
                    llh = 0:5 / 2,
                    posterior_logprob = log(c(0.5, 0.5, 2.5, 2.5, 0.5, 1.5)))
  setkeyv(lps, c("region", "event"))
  
  lse <- function(x) log(sum(exp(x)))
  regions <- data.table(location = rep(1:2, 2), region = c(1:3, 3),
                        key = "region")
  elpm <- event_logprobability_map(lps, regions)
  expect_equal(elpm[, posterior_logprob], 
               vapply(list(log(c(0.5, 0.5)),
                           log(c(0.5, 1.5)),
                           log(c(2.5, 0.5)),
                           log(c(2.5, 1.5))),
                      lse,
                      numeric(1)),
               log(1:4))
})