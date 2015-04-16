context("FastSubsetScan - general functions")


test_that("aggregate_CB_poisson: works as intended", {
  counts <- data.table(location = rep(1:2, each = 6),
                       stream = rep(1:2, 2, each = 3),
                       duration = rep(1:3, 4))
  counts[, count := 1:12]
  counts[, baseline := 1:12]
  CBim <- aggregate_CB_poisson(counts)
  expected <- c(1, 1+2, 3+3, 4, 4+5, 9+6, 7, 7+8, 15+9, 10, 10+11, 21+12)
  expect_equal(CBim[, aggregate_count], expected)
  expect_equal(CBim[, aggregate_baseline], expected)
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


test_that("aggregate_per_stream_regionasatomic: works as intended", {
  aggreg <- data.table(region = rep(3L, 8),
                       duration = rep(1:2, each = 4),
                       stream = rep(1:2, 2, each = 2),
                       location = rep(1:2, 4),
                       aggregate_count = 1:8,
                       aggregate_baseline = 1:8)

  ag <- aggregate_per_stream_regionasatomic(aggreg)
  expect_equal(ag[, aggregate_count], c(3, 7, 11, 15))
  expect_equal(ag[, aggregate_baseline], c(3, 7, 11, 15))
})

test_that("aggregate_per_stream_regionaslist: works as intended", {
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2, 4),
                    aggregate_count = 1:4,
                    aggregate_baseline = 4:1)
  locs <- 1:2
  expected_ac <- c(4, 6)
  expected_ab <- c(6, 4)
  res <- aggregate_per_stream_regionaslist(ags, locs)
  expect_equal(res[, aggregate_count], expected_ac)
  expect_equal(res[, aggregate_baseline], expected_ab)
  expect_equal(res[1, region], list(1:2))
})