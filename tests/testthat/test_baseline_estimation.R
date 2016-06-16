context("Baseline estimation")

test_that("kulldorff_baseline: calculated correctly", {
  counts <- data.table(stream = rep(1:2, each = 4),
                       location = rep(1:2, 2, each = 2),
                       time = rep(0:1, 4),
                       count = 1:8)
  setkeyv(counts, c("stream", "location", "time"))
  bs <- kulldorff_baseline(counts)
  ks <- c(c((1+2)*(1+3), 
           (1+2)*(2+4),
           (3+4)*(1+3),
           (3+4)*(2+4)) / sum(1:4),
         c((5+6)*(5+7),
           (5+6)*(6+8),
           (7+8)*(5+7),
           (7+8)*(6+8)) / sum(5:8))
  expect_equal(bs[, baseline], ks)
})