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


test_that("aggregate_CB_gaussian: works as intended", {
  counts <- data.table(location = rep(1:2, each = 4),
                       stream = rep(1:2, 2, each = 2),
                       duration = rep(1:2, 4))
  counts[, count := 1:8]
  counts[, baseline := 1:8]
  counts[, variance := 1:8]
  CBim <- aggregate_CB_gaussian(counts)
  expect_equal(CBim[, aggregate_count], c(1, 3, 3, 7, 5, 11, 7, 15))
  expect_equal(CBim[, aggregate_baseline], c(1, 3, 3, 7, 5, 11, 7, 15))
})

test_that("aggregate_CB_exponential: works as intended", {
  counts <- data.table(location = rep(1:2, each = 4),
                       stream = rep(1:2, 2, each = 2),
                       duration = rep(1:2, 4))
  counts[, count := 1:8]
  counts[, baseline := 1:8]
  CBim <- aggregate_CB_exponential(counts)
  expect_equal(CBim[, aggregate_count], rep(1:2, 4))
  expect_equal(CBim[, aggregate_baseline], rep(1:2, 4))
})


test_that("aggregate_again: works as intended", {
  aggreg <- data.table(region = rep(3L, 8),
                       duration = rep(1:2, each = 4),
                       stream = rep(1:2, 2, each = 2),
                       location = rep(1:2, 4),
                       aggregate_count = 1:8,
                       aggregate_baseline = 1:8)

  ag <- aggregate_again(aggreg)
  expect_equal(ag[, aggregate_count], c(3, 7, 11, 15))
  expect_equal(ag[, aggregate_baseline], c(3, 7, 11, 15))
})

