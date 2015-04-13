context("FSS naive Kulldorff functions")

test_that("score_minimal_stream_subset: works as intended", {
  sco <- data.table(region = rep(1:2, each = 4),
                    duration = rep(1:2, 2, each = 2),
                    stream = rep(1:2, 4),
                    score = c(-1, -1, 1, 0, 1:4))
  rr <- c(1.2, 1.3)
  expected_score <- c(1, 3, 7)
  expected_is <- list(1L, 1:2, 1:2)
  res <- score_minimal_stream_subset(sco)
  expect_equal(res[, score], expected_score)
  expect_equal(res[, included_streams], expected_is)
})