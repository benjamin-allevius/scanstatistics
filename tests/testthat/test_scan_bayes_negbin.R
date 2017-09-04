context("Bayes negbin tests")

lprob <- function(C, B, a, b) a * log(b) + lgamma(a + C) - (a+C)*log(b+B) - lgamma(a)
prob <- function(C, B, a, b) exp(lprob(C, B, a, b))

beg_manual <- function(counts, baselines, zones, outbreak_prior, alpha_null, beta_null,
                       alpha_alt, beta_alt) {
 
  
  win_post <- data.frame(zone = rep(NA, length(zones) * nrow(counts)))
  win_post$duration <- NA
  win_post$posterior <- NA
  idx <- 1
  for (j in 1:nrow(counts)) {
    for (i in seq_along(zones)) {
      z <- zones[[i]]
      win_post[idx, 1] <- i
      win_post[idx, 2] <- j
      win_post[idx, 3] <- prob(sum(counts[1:j, z]), 
                               sum(baselines[1:j, z]), 
                               alpha_alt, beta_alt) * 
                          prob(sum(counts) - sum(counts[1:j, z]), 
                               sum(baselines) - sum(baselines[1:j, z]), 
                               alpha_null, beta_null) * 
                          outbreak_prior  / (length(zones) * nrow(counts))
      idx <- idx + 1
    }
  }
  
  null_post <- prob(sum(counts), sum(baselines), alpha_null, beta_null) * (1 - outbreak_prior)
  
  p_d <- sum(win_post$posterior) + null_post
  win_post$posterior <- win_post$posterior / p_d
  
  win_post$bayes_factor <- win_post$posterior / null_post
  
  list(null_prior = 1 - outbreak_prior,
       null_post = null_post / p_d,
       outbreak_prior = outbreak_prior,
       outbreak_post = sum(win_post$posterior),
       marginal_data_prob = p_d,
       window_prior = outbreak_prior  / (length(zones) * nrow(counts)),
       window_post = win_post)
  
}

test_that("scan_bayes_negbin", {


  # Single timepoint -----------------------------------------------------------
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 2), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))
  
  outbreak_prior <- 0.1
  
  
  expected <- beg_manual(in1$counts,
                         in1$baselines,
                         in1$zones,
                         outbreak_prior,
                         2, 2,
                         4, 2)

  actual <- scan_bayes_negbin_cpp(in1$counts,
                                   in1$baselines,
                                   in1$zones_flat - 1,
                                   in1$zone_lengths,
                                   outbreak_prior,
                                   alpha_null = 2,
                                   beta_null = 2,
                                   alpha_alt = 4,
                                   beta_alt = 2)

  # 3 timepoints ---------------------------------------------------------------
  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    baselines = matrix(c(0.5, 2,
                         0.5, 2,
                         0.5, 2), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))
  
  expected <- beg_manual(in2$counts,
                         in2$baselines,
                         in2$zones,
                         outbreak_prior,
                         2, 2,
                         4, 2)

  actual <- scan_bayes_negbin_cpp(in2$counts,
                                  in2$baselines,
                                  in2$zones_flat - 1,
                                  in2$zone_lengths,
                                  0.1,
                                  alpha_null = 2,
                                  beta_null = 2,
                                  alpha_alt = 4,
                                  beta_alt = 2)
})
