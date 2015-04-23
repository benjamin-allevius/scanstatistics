context("FastSubsetScan - naive Kulldorff method")

test_that("naive_kulldorff: works as intended", {
  cts <- table_creator(list(location = 1:2, duration = 1:2, stream = 1:2))
  cts[, count := 6]
  cts[, baseline := 5]
  setkeyv(cts, c("location", "duration", "stream"))
  regions <- sets::set(sets::as.set(1L),sets::as.set(2L),sets::as.set(1:2))
  # As counts are always higher than baselines, region 3 (locations 1 and 2)
  # with streams 1 and 2, and duration 2, should be the top-scoring combination
  res <- naive_kulldorff(cts, "poisson", regions)
  setkey(res, "score")
  expect_equal(res[.N, region], 3L)
  expect_equal(res[.N, duration], 2L)
  expect_equal(res[.N, included_streams], list(1:2))
})

test_that("naive_kulldorff: works as intended", {
  cts <- table_creator(list(location = 1:2, duration = 1:2, stream = 1:2))
  cts[, count := 5]
  cts[, baseline := 6]
  setkeyv(cts, c("location", "duration", "stream"))
  regions <- sets::set(sets::as.set(1L),sets::as.set(2L),sets::as.set(1:2))
  # As baselines are always higher than counts, we should get an empty table
  # as results
  res <- naive_kulldorff(cts, "poisson", regions)
  expect_equal(nrow(res), 0L)
})