context("FastSubsetScan - fast Kulldorff method")

# find_maximizing_subsets -----------------------------------------------------------------

test_that("find_maximizing_subsets: works as intended", {
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2L, 4),
                    aggregate_count = 1:4,
                    aggregate_baseline = 4:1)
  score_fun <- score_fun_EBP
  cond_score_fun <- conditional_score_fun_EBP
  expected_score <- score_fun(3, 2) + score_fun(4, 1)
  expected_streams <- list(1:2)
  expected_reg <- list(2L)
  res <- find_maximizing_subsets(aggregates = ags, 
                                 score_fun = score_fun,
                                 cond_score_fun = cond_score_fun)
  expect_equal(res[, score], expected_score)
  expect_equal(res[, included_streams], expected_streams)
  expect_equal(res[, region], expected_reg)
})

test_that("find_maximizing_subsets: works as intended with 'bad' input", {
  # counts are equal to baselines, so no subset of of locations or streams
  # are optimal
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2L, 4),
                    aggregate_count = 4:1,
                    aggregate_baseline = 4:1)
  score_fun <- score_fun_EBP
  cond_score_fun <- conditional_score_fun_EBP
  res <- find_maximizing_subsets(aggregates = ags, 
                                 score_fun = score_fun,
                                 cond_score_fun = cond_score_fun)
  expect_equal(res[, score], as.numeric(NA))
  expect_equal(res[, included_streams], list(NULL))
  expect_equal(res[, region], list(NULL))
})

# optimal_stream_subset --------------------------------------------------------

test_that("optimal_stream_subset: works as intended", {
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2L, 4),
                    aggregate_count = 1:4,
                    aggregate_baseline = 4:1)
  score_fun <- score_fun_EBP
  expected_score <- score_fun(3, 2) + score_fun(4, 1)
  expected_streams <- list(1:2)
  res <- optimal_stream_subset(ags, 2L, score_fun)
  expect_equal(res[, score], expected_score)
  expect_equal(res[, included_streams], expected_streams)
})

test_that("optimal_stream_subset: works as intended with empty output", {
  # counts are equal to baselines, so no subset of streams is optimal
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2L, 4),
                    aggregate_count = 4:1,
                    aggregate_baseline = 4:1)
  score_fun <- score_fun_EBP
  res <- optimal_stream_subset(ags, 2L, score_fun)
  expect_equal(nrow(res), 0)
})

# find_maximizing_region -------------------------------------------------------

test_that("find_maximizing_region: works as intended", {
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2L, 4),
                    aggregate_count = 1:4,
                    aggregate_baseline = 4:1)
  cond_score_fun <- conditional_score_fun_EBP
  expected <- c(cond_score_fun(1, 4, 1.2) + 
                  cond_score_fun(2, 3, 1.3),
                cond_score_fun(3, 2, 1.2) +
                  cond_score_fun(4, 1, 1.3))
  res <- find_maximizing_region(aggregates = ags, 
                                cond_score_fun = cond_score_fun, tol = 0.01)
  expect_equal(res, 2L)
})


test_that("fast_kulldorff_priority: works as intended", {
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2L, 4),
                    aggregate_count = 1:4,
                    aggregate_baseline = 4:1)
  rr <- c(1.2, 1.3)
  cond_score_fun <- conditional_score_fun_EBP
  expected <- c(cond_score_fun(1, 4, 1.2) + 
                  cond_score_fun(2, 3, 1.3),
                cond_score_fun(3, 2, 1.2) +
                  cond_score_fun(4, 1, 1.3))
  res <- fast_kulldorff_priority(ags, rr, cond_score_fun)
  expect_equal(res[, priority], expected)
  expect_equal(res[1, included_streams], list(1:2))
})

test_that("fast_kulldorff_maxregion: works as intended", {
  pri <- data.table(location = 1:5,
                    duration = rep(2L, 5),
                    included_streams = list(1:2, 1:2, 1:2, 1:2, 1:2),
                    priority = c(-0.5, 0.5, 2, 0, 3.4))
  expected <- 0.5 + 2 + 3.4
  res <- fast_kulldorff_maxregion(pri)
  expect_equal(res[, score], expected)
  expect_equal(res[1, region], list(c(2L, 3L, 5L)))
})

test_that("relative_risk_mle: works as intended", {
  # Recall that for Fast Kulldorff, we don't have a column for region,
  # instead we work with the locations directly, and pick out those that
  # belong to the region we want.
  ags <- data.table(location = rep(1:3, each = 3),
                    stream = rep(1:3, 3),
                    duration = rep(2L, 9), 
                    key = "stream")
  # Aggregates are the C_i^m(W), B_i^m(W)
  ags[, aggregate_count := 1:9]
  ags[, aggregate_baseline := 9:1]
  cs <- c(2+3, 5+6, 8+9)
  bs <- c(8+7, 5+4, 2+1)
  rrs <- cs / bs
  expected <- c(1, rrs[2:3]) 
  res <- relative_risk_mle(ags, c(2L, 3L)) # pick region = locations 2, 3
  expect_equal(unname(res), expected)
  expect_true(all(as.integer(names(res)) == 1:3))
})