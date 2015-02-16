#library(data.table)
context("Event probability map")


# locations 1, 2, 3, 4
# regions (1), (2), (3), (4), (1,2), (2,4), (1,3,4)
nevs <- 2
repp <- data.table(region = rep(c(1:5, c(5, 6, 6, 7, 7, 7)), nevs),
                   event = rep(1:nevs, each = 11),
                   location = rep(c(1:4, c(1, 2, 2, 4, 1, 3, 4)), nevs),
                   probability = 0)
setkey(repp, location)
repp[, "probability"] <- c(rep(1:4, each = 3) / 3,
                         rep(5:6, each = 2) / 2,
                         rep(7:8, each = 3) / 3)

test_that("output data.table has correct names", {
    pm <- event_probability_map(repp)
    expect_equal(names(pm), c("event", "location", "probability"))
})


test_that("probabilities sum correctly", {
    pm <- event_probability_map(repp)
    expect_equal(pm[, probability], 1:8)
})