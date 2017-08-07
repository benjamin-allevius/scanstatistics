context("Space-time permutation statistic")

expand_matrix <- function(A) {
  res <- matrix(NA, nrow = sum(A), ncol = 2)
  index <- 1
  for (j in 1:ncol(A)) {
    for (i in 1:nrow(A)) {
      n <- A[i,j]
      for (k in seq_len(n)) {
        res[index, 1] <- i
        res[index, 2] <- j
        index <- index + 1
      }
    }
  }
  return(res)
}

contract_matrix <- function(A, nr, nc) {
  res <- matrix(0L, nr, nc)
  for (i in 1:nrow(A)) {
    res[A[i, 1], A[i, 2]] <- res[A[i, 1], A[i, 2]] + 1L
  }
  return(res)
}


permute_table <- function(A) {
  a <- expand_matrix(A)
  for (i in nrow(a):1) {
    j <- floor(runif(1, 1, nrow(a)))
    tmp <- a[i, 1]
    a[i, 1] <- a[j, 1]
    a[j, 1] <- tmp
  }
  return(contract_matrix(a, nrow(A), ncol(A)))
}

test_that("scan_pb_perm_cpp", {

  # 3 timepoints
  in2 <- list(
    counts = matrix(c(1, 0,
                      2, 1,
                      0, 20), nrow = 3, byrow = TRUE),
    zones = list(1L, 2L, 1:2))
  in2$N <- sum(in2$counts)
  in2$baselines <- outer(rowSums(in2$counts), colSums(in2$counts)) / in2$N
  in2$zones_flat =  unlist(in2$zones)
  in2$zone_lengths = unlist(lapply(in2$zones, length))
  

  perm1 <- scan_pb_perm_cpp(in2$counts,
                              in2$baselines,
                              in2$zones_flat - 1,
                              in2$zone_lengths,
                              store_everything = FALSE,
                              num_mcsim = 0)$observed
  pois1 <- scan_pb_poisson_cpp(in2$counts,
                            in2$baselines,
                            in2$zones_flat - 1,
                            in2$zone_lengths,
                            store_everything = FALSE,
                            num_mcsim = 0)$observed
  
  expect_equal(perm1, pois1)
  
  set.seed(1)
  in2$counts <- permute_matrix(in2$counts)
  
  # pois2 <- scan_pb_poisson_cpp(in2$counts,
  #                              in2$baselines,
  #                              in2$zones_flat - 1,
  #                              in2$zone_lengths,
  #                              store_everything = FALSE,
  #                              num_mcsim = 0)$observed
  # set.seed(1)
  # perm2 <- scan_pb_perm_cpp(in2$counts,
  #                           in2$baselines,
  #                           in2$zones_flat - 1,
  #                           in2$zone_lengths,
  #                           store_everything = FALSE,
  #                           num_mcsim = 1)$simulated[1, ]
  # 
  # expect_equal(perm2, pois2)
})
