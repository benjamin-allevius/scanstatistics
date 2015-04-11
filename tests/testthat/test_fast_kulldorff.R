context("FSS fast Kulldorff functions")


test_that("fast_kulldorff_priority: works as intended", {
  ags <- data.table(location = rep(1:2, each = 2),
                    stream = rep(1:2, 2),
                    duration = rep(2, 4),
                    aggregate_count = 1:4,
                    aggregate_baseline = 4:1)
  rr <- c(1.2, 1.3)
  expected <- c(fk_priority_term_poisson(1, 4, 1.2) + 
                  fk_priority_term_poisson(2, 3, 1.3),
                fk_priority_term_poisson(3, 2, 1.2) +
                  fk_priority_term_poisson(4, 1, 1.3))
  res <- fast_kulldorff_priority(ags, rr, fk_priority_term_poisson)
  expect_equal(res[, priority], expected)
  expect_equal(res[1, included_streams], list(1:2))
})

test_that("fast_kulldorff_maxregion: works as intended", {
  pri <- data.table(location = 1:5,
                    duration = rep(2, 5),
                    priority = c(-0.5, 0.5, 2, 0, 3.4),
                    included_streams = list(1:2, 1:2, 1:2, 1:2, 1:2))
  expected <- 0.5 + 2 + 3.4
  res <- fast_kulldorff_maxregion(pri)
  expect_equal(res[, score], expected)
  expect_equal(res[1, included_locations], list(c(2L, 3L, 5L)))
})