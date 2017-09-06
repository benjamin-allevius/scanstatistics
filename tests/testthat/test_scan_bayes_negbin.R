context("Bayes negbin tests")

lprob <- function(C, B, a, b) {
  a * log(b) + lgamma(a + C) - (a+C)*log(b+B) - lgamma(a) - lgamma(1 + C)
}

prob <- function(C, B, a, b) exp(lprob(C, B, a, b))

win_score <- function(C_in, B_in, a_in, b_in, C_out, B_out, a_out, b_out, pri) {
  prob(C_in, B_in, a_in, b_in) * prob(C_out, B_out, a_out, b_out) * pri
}

beg_manual <- function(counts, baselines, zones, 
                       alt_prior, 
                       alpha_null, beta_null,
                       alpha_alt, beta_alt,
                       inc_vals, inc_probs) {
 
  
  win_post <- data.frame(zone = rep(NA, length(zones) * nrow(counts)))
  win_post$duration <- NA
  win_post$posterior <- 0
  
  win_cpost <- matrix(0, nrow(win_post), length(inc_vals))
  
  null_cpost <- rep(0, length(inc_vals))
  alt_cpost <- rep(0, length(inc_vals))
  m_post <- rep(0, length(inc_vals))
  data_cpost <- rep(0, length(inc_vals))
  
  
  for (k in 1:length(inc_vals)) {
    idx <- 1
    for (j in 1:nrow(counts)) {
      for (i in seq_along(zones)) {
        z <- zones[[i]]
        win_post[idx, "zone"] <- i
        win_post[idx, "duration"] <- j
        win_cpost[idx, k] <- win_score(
                          sum(counts[1:j, z]),
                          sum(baselines[1:j, z]),
                          inc_vals[k] * alpha_alt, beta_alt,
                          sum(counts) - sum(counts[1:j, z]),
                          sum(baselines) - sum(baselines[1:j, z]),
                          alpha_null, beta_null,
                          alt_prior / (length(zones) * nrow(counts)))
        idx <- idx + 1
      }
    }
    np <- prob(sum(counts), sum(baselines), 
               alpha_null, beta_null) * (1 - alt_prior)
    data_cpost[k] <- sum(win_cpost[, k]) + np
    win_cpost[, k] <- win_cpost[, k] / data_cpost[k]
    null_cpost[k] = np / data_cpost[k]
    alt_cpost[k] = 1 - null_cpost[k]
    m_post[k] <- data_cpost[k] * inc_probs[k]
  }
  
  pd <- sum(data_cpost * inc_probs)
  
  null_post <- sum(null_cpost * m_post) / pd
  alt_post <- 1 - null_post
  
  m_post <- data_cpost * inc_probs / pd
  
  for (k in 1:length(inc_vals)) {
    win_post$posterior <- win_post$posterior + win_cpost[, k] * m_post[k]
  }
  win_post$log_posterior <- log(win_post$posterior)
  
  list(priors = 
         list(null_prior = 1 - alt_prior,
              alt_prior = alt_prior,
              window_prior = alt_prior  / (length(zones) * nrow(counts)),
              m_prior = inc_probs),
       posteriors =
         list(null_posterior = null_post,
              alt_posterior = alt_post,
              inc_posterior = m_post,
              window_posteriors = win_post),
       marginal_data_prob = pd)
  
}

test_that("scan_bayes_negbin: 1 timepoint", {


  # Single m -----------------------------------------------------------
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 2), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))
  
  expected <- beg_manual(in1$counts,
                         in1$baselines,
                         in1$zones,
                         0.1,
                         2, 2,
                         4, 2,
                         1, 1)

  actual <- scan_bayes_negbin_cpp(in1$counts,
                                  in1$baselines,
                                  in1$zones_flat - 1,
                                  in1$zone_lengths,
                                  0.1,
                                  alpha_null = 2,
                                  beta_null = 2,
                                  alpha_alt = 4,
                                  beta_alt = 2,
                                  inc_values = 1,
                                  inc_probs = 1)
  
  expect_equal(actual$posteriors$null_posterior,
               expected$posteriors$null_posterior)
  expect_equal(actual$posteriors$alt_posterior,
               expected$posteriors$alt_posterior)
  expect_equal(as.vector(actual$posteriors$inc_posterior$inc_posterior),
               expected$posteriors$inc_posterior)
  expect_equal(actual$posteriors$window_posteriors$log_posterior,
               expected$posteriors$window_posteriors$log_posterior)
  
  
  # Multiple values of m -------------------------------------------------------
  inc_vals <- c(1, 2)
  inc_probs = c(0.5, 0.5)
  
  expected <- beg_manual(in1$counts,
                         in1$baselines,
                         in1$zones,
                         0.1,
                         2, 2,
                         4, 2,
                         inc_vals, inc_probs)
  
  actual <- scan_bayes_negbin_cpp(in1$counts,
                                  in1$baselines,
                                  in1$zones_flat - 1,
                                  in1$zone_lengths,
                                  0.1,
                                  alpha_null = 2,
                                  beta_null = 2,
                                  alpha_alt = 4,
                                  beta_alt = 2,
                                  inc_values = inc_vals,
                                  inc_probs = inc_probs)
  
  expect_equal(actual$posteriors$null_posterior,
               expected$posteriors$null_posterior)
  expect_equal(actual$posteriors$alt_posterior,
               expected$posteriors$alt_posterior)
  expect_equal(as.vector(actual$posteriors$inc_posterior$inc_posterior),
               expected$posteriors$inc_posterior)
  expect_equal(actual$posteriors$window_posteriors$log_posterior,
               expected$posteriors$window_posteriors$log_posterior)
})


test_that("scan_bayes_negbin: 3 timepoints", {
  
  # 1 value of m ---------------------------------------------------------------
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
                         0.1,
                         2, 2,
                         4, 2,
                         1, 1)
  
   actual <- scan_bayes_negbin_cpp(in2$counts,
                                  in2$baselines,
                                  in2$zones_flat - 1,
                                  in2$zone_lengths,
                                  0.1,
                                  alpha_null = 2,
                                  beta_null = 2,
                                  alpha_alt = 4,
                                  beta_alt = 2,
                                  inc_values = 1,
                                  inc_probs = 1)
   
   expect_equal(actual$posteriors$null_posterior,
                expected$posteriors$null_posterior)
   expect_equal(actual$posteriors$alt_posterior,
                expected$posteriors$alt_posterior)
   expect_equal(as.vector(actual$posteriors$inc_posterior$inc_posterior),
                expected$posteriors$inc_posterior)
   expect_equal(actual$posteriors$window_posteriors$log_posterior,
                expected$posteriors$window_posteriors$log_posterior)
  
  # Multiple values of m -------------------------------------------------------
   inc_vals <- seq(1, 3, 0.1)
   inc_probs <- rep(1 / length(inc_vals), length(inc_vals))
   
   expected <- beg_manual(in2$counts,
                          in2$baselines,
                          in2$zones,
                          0.1,
                          2, 2,
                          4, 2,
                          inc_vals, inc_probs)
   
   actual <- scan_bayes_negbin_cpp(in2$counts,
                                   in2$baselines,
                                   in2$zones_flat - 1,
                                   in2$zone_lengths,
                                   0.1,
                                   alpha_null = 2,
                                   beta_null = 2,
                                   alpha_alt = 4,
                                   beta_alt = 2,
                                   inc_values = inc_vals,
                                   inc_probs = inc_probs)
   
   expect_equal(actual$posteriors$null_posterior,
                expected$posteriors$null_posterior)
   expect_equal(actual$posteriors$alt_posterior,
                expected$posteriors$alt_posterior)
   expect_equal(as.vector(actual$posteriors$inc_posterior$inc_posterior),
                expected$posteriors$inc_posterior)
   expect_equal(actual$posteriors$window_posteriors$log_posterior,
                expected$posteriors$window_posteriors$log_posterior)
})
