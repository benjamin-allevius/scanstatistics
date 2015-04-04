context("FSS score functions")


test_that("aggregate_CB_poisson: works as intended", {
  counts <- data.table(location = rep(1:2, each = 4),
                       stream = rep(1:2, 2, each = 2),
                       duration = rep(1:2, 4))
  counts[, count := 1:8]
  counts[, baseline := 1:8]
  CBim <- aggregate_CB_poisson(counts)
  expect_equal(CBim[, aggregate_count], c(1, 3, 3, 7, 5, 11, 7, 15))
  expect_equal(CBim[, aggregate_baseline], c(1, 3, 3, 7, 5, 11, 7, 15))
})






