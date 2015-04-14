context("FastSubsetScan - fast Kulldorff functions")


test_that("fast_kulldorff_priority: works as intended", {
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2, 4),
                    aggregate_count = 1:4,
                    aggregate_baseline = 4:1)
  rr <- c(1.2, 1.3)
  expected <- c(priority_term_poisson(1, 4, 1.2) + 
                  priority_term_poisson(2, 3, 1.3),
                priority_term_poisson(3, 2, 1.2) +
                  priority_term_poisson(4, 1, 1.3))
  res <- fast_kulldorff_priority(ags, rr, priority_term_poisson)
  expect_equal(res[, priority], expected)
  expect_equal(res[1, included_streams], list(1:2))
})

test_that("fast_kulldorff_maxregion: works as intended", {
  pri <- data.table(location = 1:5,
                    duration = rep(2, 5),
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
                    duration = rep(2, 9), 
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