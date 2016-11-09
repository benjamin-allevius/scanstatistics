context("FSS common functions")

# test_that("aggregate_locations", {
#   table <- data.table(expand.grid(location = 1:2,
#                                   duration = 1:2,
#                                   stream = 1:3))
#   table[, mu := location /  duration + stream]
#   table[, count := rpois(.N, mu)]
#   
#   stream_subsets <- list(1, 2, 3, c(1, 2), c(1, 3))
#   
#   
#   expect_equal(actual, expected)
# })