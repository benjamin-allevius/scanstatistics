context("FastSubsetScan - general functions")

# score_minimal_stream_subset---------------------------------------------------

test_that("score_minimal_stream_subset: works as intended for atomic region", {
  sco <- data.table(region = rep(1:2, each = 4),
                    duration = rep(1:2, 2, each = 2),
                    stream = rep(1:2, 4),
                    score = c(-1, -1, 1, 0, 1:4))
  rr <- c(1.2, 1.3)
  expected_score <- c(1, 3, 7)
  expected_is <- list(1L, 1:2, 1:2)
  res <- score_minimal_stream_subset(sco, region_as_list = FALSE)
  expect_equal(res[, score], expected_score)
  expect_equal(res[, included_streams], expected_is)
})

test_that("score_minimal_stream_subset: works as intended for list region", {
  sco <- data.table(region = rep(list(1:2), 4),
                    duration = rep(1L, 4),
                    stream = 1:4,
                    score = c(-1, 1, 0, 1))
  expected_score <- 2
  expected_region <- list(1:2)
  expected_streams <- list(c(2L, 4L))
  res <- score_minimal_stream_subset(sco, region_as_list = TRUE)
  expect_equal(res[, score], expected_score)
  expect_equal(res[, region], expected_region)
  expect_equal(res[, included_streams], expected_streams)
})

test_that("score_minimal_stream_subset_regionasatomic: works as intended", {
  sco <- data.table(region = rep(1:2, each = 4),
                    duration = rep(1:2, 2, each = 2),
                    stream = rep(1:2, 4),
                    score = c(-1, -1, 1, 0, 1:4))
  rr <- c(1.2, 1.3)
  expected_score <- c(1, 3, 7)
  expected_is <- list(1L, 1:2, 1:2)
  res <- score_minimal_stream_subset_regionasatomic(sco)
  expect_equal(res[, score], expected_score)
  expect_equal(res[, included_streams], expected_is)
})

test_that("score_minimal_stream_subset_regionaslist: works as intended", {
  sco <- data.table(region = rep(list(1:2), 4),
                    duration = rep(1L, 4),
                    stream = 1:4,
                    score = c(-1, 1, 0, 1))
  expected_score <- 2
  expected_region <- list(1:2)
  expected_streams <- list(c(2L, 4L))
  res <- score_minimal_stream_subset_regionaslist(sco)
  expect_equal(res[, score], expected_score)
  expect_equal(res[, region], expected_region)
  expect_equal(res[, included_streams], expected_streams)
})

# aggregate_per_stream ---------------------------------------------------------

test_that("aggregate_per_stream: works as intended  for atomic region", {
  aggreg <- data.table(region = rep(3L, 8),
                       duration = rep(1:2, each = 4),
                       stream = rep(1:2, 2, each = 2),
                       location = rep(1:2, 4),
                       aggregate_count = 1:8,
                       aggregate_baseline = 1:8)
  
  ag <- aggregate_per_stream(aggreg, region_as_list = FALSE)
  expect_equal(ag[, aggregate_count], c(3, 7, 11, 15))
  expect_equal(ag[, aggregate_baseline], c(3, 7, 11, 15))
})

test_that("aggregate_per_stream: works as intended  for list region", {
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2, 4),
                    aggregate_count = 1:4,
                    aggregate_baseline = 4:1)
  locs <- 1:2
  expected_ac <- c(4, 6)
  expected_ab <- c(6, 4)
  res <- aggregate_per_stream(ags, locs, region_as_list = TRUE)
  expect_equal(res[, aggregate_count], expected_ac)
  expect_equal(res[, aggregate_baseline], expected_ab)
  expect_equal(res[1, region], list(1:2))
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

# aggregate_CB -----------------------------------------------------------------

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