context("MBSS main functions")

test_that("probability_map: works for valid input", {
  
  llrs <- table_creator(list(event = 1:2,
                             stream = 1:2,
                             location = 1:2, 
                             time = 0:2))
  set.seed(123)
  cts <- rpois(nrow(llrs), lambda = 8)
  null_densities <- dpois(cts, lambda = 8, log = TRUE)
  alt_densities <- dpois(cts, lambda = 10, log = TRUE)
  null_llh <- sum(null_densities)
  
  llrs[, llr := alt_densities - null_densities]
  
  null_prior <- 0.9
  event_priors <- c(0.03, 0.07)
  duration_condpriors <- matrix(c(1:3 / 6, 3:1 / 6), ncol = 2)
  
  regions <- sets::set(sets::as.set(1L), sets::as.set(2L), sets::as.set(1:2))
  
  m <- MBSS(llrs, 
            null_llh, 
            regions, 
            null_prior, 
            event_priors, 
            duration_condpriors)
})