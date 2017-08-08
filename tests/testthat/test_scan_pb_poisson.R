context("PB Poisson statistic tests")

pb_rrin <- function(C,B) sum(C) / sum(B)
pb_rrout <- function(C,B,N) {
  C <- sum(C); B <- sum(B) 
  ifelse(B < N, (N - C) / (N-B), 1)
}

pb_score <- function(C, B, N) {
  C <- sum(C)
  B <- sum(B)
  risk_in <- C / B
  risk_out <- ifelse(N > B, (N - C) / (N - B), 1)
  term2 <- ifelse(N > C, (N-C)*log(risk_out), 0)
  ifelse(C > B, C * log(risk_in) + term2, 0)
}



test_that("scan_pb_poisson_cpp", {
  # Single timepoint
  in1 <- list(
    counts = matrix(c(1, 0), nrow = 1),
    baselines = matrix(c(0.5, 0.5), nrow = 1),
    zones = list(1L, 2L, 1:2))
  in1$N = sum(in1$counts)
  in1$zones_flat =  unlist(in1$zones)
  in1$zone_lengths = unlist(lapply(in1$zones, length))
  
  ex1_score <- rep(NA, length(in1$zones) * nrow(in1$counts))
  ex1_relrisk_in <- ex1_score
  ex1_relrisk_out <- ex1_score
  idx <- 1
  for (j in 1:nrow(in1$counts)) {
    for (z in in1$zones) {
      ex1_score[idx] <- pb_score(in1$counts[1:j, z], in1$baselines[1:j, z], in1$N)
      ex1_relrisk_in[idx] <- pb_rrin(in1$counts[1:j, z], in1$baselines[1:j, z])
      ex1_relrisk_out[idx] <- pb_rrout(in1$counts[1:j, z], in1$baselines[1:j, z], in1$N)
      idx <- idx + 1
    }
  }
  

  actual1 <- scan_pb_poisson_cpp(in1$counts,
                                 in1$baselines,
                                 in1$zones_flat - 1,
                                 in1$zone_lengths,
                                 store_everything = TRUE,
                                 num_mcsim = 0)$observed
  actual1b <- scan_pb_poisson_cpp(in1$counts,
                                  in1$baselines,
                                 in1$zones_flat - 1,
                                 in1$zone_lengths,
                                 store_everything = FALSE,
                                 num_mcsim = 0)$observed
  expect_equal(actual1$score, ex1_score)
  expect_equal(actual1$relrisk_in, ex1_relrisk_in)
  expect_equal(actual1$relrisk_out, ex1_relrisk_out)
  expect_equal(c(actual1[which.max(actual1$score), ]), c(actual1b))

  # 3 timepoints
  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    baselines = matrix(c(2, 1.0,
                         4, 5,
                         6, 6), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))
  in2$N <- sum(in2$counts)
  
  ex2_score <- rep(NA, length(in2$zones) * nrow(in2$counts))
  ex2_relrisk_in <- ex2_score
  ex2_relrisk_out <- ex2_score
  idx <- 1
  for (j in 1:nrow(in2$counts)) {
    for (z in in2$zones) {
      ex2_score[idx] <- pb_score(in2$counts[1:j, z], in2$baselines[1:j, z], in2$N)
      ex2_relrisk_in[idx] <- pb_rrin(in2$counts[1:j, z], in2$baselines[1:j, z])
      ex2_relrisk_out[idx] <- pb_rrout(in2$counts[1:j, z], in2$baselines[1:j, z], in2$N)
      idx <- idx + 1
    }
  }
  
  actual2 <- scan_pb_poisson_cpp(in2$counts,
                                 in2$baselines,
                                 in2$zones_flat - 1,
                                 in2$zone_lengths,
                                 store_everything = TRUE,
                                 num_mcsim = 0)$observed
  actual2b <- scan_pb_poisson_cpp(in2$counts,
                                  in2$baselines,
                                 in2$zones_flat - 1,
                                 in2$zone_lengths,
                                 store_everything = FALSE,
                                 num_mcsim = 0)$observed

  expect_equal(actual2$score, ex2_score)
  expect_equal(actual2$relrisk_in, ex2_relrisk_in)
  expect_equal(actual2$relrisk_out, ex2_relrisk_out)
  expect_equal(c(actual2[which.max(actual2$score), ]), c(actual2b))
})

test_that("scan_pb_poisson_cpp mcsim", {

  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    baselines = matrix(c(2, 1.0,
                         4, 5,
                         6, 6), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))
  in2$N <- sum(in2$counts)

  # Generate counts with rnbinom with R, run scan_pb_poisson_cpp
  set.seed(1)
  num_sims <- 2
  obs_res <- NULL
  for (i in 1:num_sims) {
    in2$counts <- matrix(rmultinom(1, in2$N, as.vector(in2$baselines)),
                       nrow(in2$baselines), ncol(in2$baselines))
    res <- scan_pb_poisson_cpp(in2$counts,
                               in2$baselines,
                               in2$zones_flat - 1,
                               in2$zone_lengths,
                               store_everything = FALSE,
                               num_mcsim = 0)$observed
    obs_res <- rbind(obs_res, res)
  }

  # Run equal number of mcsim
  set.seed(1)
  sim_res <- scan_pb_poisson_cpp(in2$counts,
                                in2$baselines,
                                in2$zones_flat - 1,
                                in2$zone_lengths,
                                store_everything = FALSE,
                                num_mcsim = 2)$simulated

  expect_equal(obs_res, sim_res)
})
