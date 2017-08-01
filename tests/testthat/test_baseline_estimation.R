context("Baseline estimation")

test_that("estimate_baselines: bad input", {
  
  # Fewer locations implied by counts than population:
  expect_error(estimate_baselines(matrix(0, 3, 2), 1:4),
               paste0("The number of locations implied must be the same in the",
                      " counts and the population arguments."))
  
  # Fewer locations implied by population than counts:
  expect_error(estimate_baselines(matrix(0, 3, 4), 1:2),
               paste0("The number of locations implied must be the same in the",
                      " counts and the population arguments."))
  
  # Counts imply spatial analysis, population space-time:
  expect_error(estimate_baselines(1:5, matrix(1, 2, 5)),
               paste0("If counts is a vector, population should be too."))
  
})

test_that("estimate_baselines: works", {
  # Spatial analysis:
  expect_equal(estimate_baselines(1:5, 11:15),
               matrix(11:15 / sum(11:15) * sum(1:5), nrow = 1))
  # Space-time analysis with constant population:
  expect_equal(estimate_baselines(matrix(1:8, 4, 2), c(10, 20)),
               matrix(rep(c(10, 20), each = 4), 4, 2) / 30 * sum(1:8) / 4)
  # Space-time analysis with population varying over time:
  co <- matrix(1:4, 2, 2)
  pop <- matrix(c(1, 4, 3, 6), 2, 2)
  actual <- estimate_baselines(co, pop)
  expected <- matrix(c(1/4, 4/10, 3/4, 6/10), 2, 2) * sum(1:4) / 2
  expect_equal(actual, expected)
})
