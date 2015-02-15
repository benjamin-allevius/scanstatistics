#library(data.table)
context("Probability map")

test_that("output data.table has correct names", {
    epm <- data.table(event = rep(1:2, each = 5),
                      location = rep(1:5, 2),
                      probability = 1:10)
    pm <- probability_map(epm)
    expect_equal(names(pm), c("location", "probability"))
})


test_that("probabilities sum correctly", {
    epm <- data.table(event = rep(1:2, each = 5),
                      location = rep(1:5, 2),
                      probability = 1:10)
    pm <- probability_map(epm)
    expect_equal(pm[, probability], seq(from = 7, to = 15, by = 2))
})


test_that("handles NaN probability", {
    epm <- data.table(event = rep(1:2, each = 5),
                      location = rep(1:5, 2),
                      probability = c(NaN, 2:9, NaN))
    pm <- probability_map(epm)
    expect_equal(pm[, probability], seq(from = 7, to = 15, by = 2))
})


test_that("handles NA probability", {
    epm <- data.table(event = rep(1:2, each = 5),
                      location = rep(1:5, 2),
                      probability = c(NA, 2:9, NA))
    pm <- probability_map(epm)
    expect_equal(pm[, probability], seq(from = 7, to = 15, by = 2))
})