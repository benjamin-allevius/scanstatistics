context("MBSS likelihood ratio functions")


# Copies of functions in utility_functions.R, in case these change

# table_creator
tc <- function(col_list, key = NULL) {
  data.table(do.call(expand.grid, col_list), key = key)
}

# region_table_creator
rtc <- function(regions, key = NULL) {
  
  region_names <- names(regions)
  if (is.null(region_names)) {
    region_names <- seq_along(regions)
  }
  
  data.table(location = unlist(regions, use.names = FALSE),
             region = rep(region_names, 
                          vapply(regions, length, integer(1))),
             key = key)
}

# region_joiner
rj <- function(locations_etc, regions) {
  merge(x = rtc(regions, key = "location"), 
        y = locations_etc,
        by = "location", allow.cartesian = TRUE)
}

# Log-likelihood ratios: event of type k and severity l vs no event.
# For all locations and times, log-likelihood ratios are specified in terms 
# of one part that only depends on the data stream m, 
# and event type and severity through impact factor x;
# and another part which also depends on the location i.


test_that("sum LLR over stream is correct", {
  kc <- c("stream", "time", "severity", "event")
  
  # Part of LLR dependent on region
  DT <- tc(list(stream = 1:2,
                time = 0:1,
                severity = 1,
                event = 1, 
                region = 1:3))
  setkeyv(DT, kc)
  DT[, slg := rep(1:6, 2) / 2]
  
  # Part of LLR independent of region
  DT2 <- tc(list(stream = 1:2,
                 time = 0:1,
                 severity = 1,
                 event = 1),
            key = kc)
  DT2[, lf := c(-2, -1, 2, 1)]
  
  stream_term <- sum_over_streams(DT, DT2)
  expect_equal(stream_term[, lq], 1:6)
})

# sum_locations_in_region
test_that("sum over locations in region", {
  
  DT <- tc(list(stream = 1:2,
                time = 0:1,
                severity = 1,
                event = 1, 
                location = 1:2))
  DT <- rj(DT, list(1, 2, 1:2))
  
  kc <- c("stream", "time", "severity", "event")
  setkeyv(DT, c("region", kc))
  
  DT[, lg := c(1:8, rep(c(10, 100, 1000, 10000) / 2, each = 2))]
  
  region_sums <- sum_locations_in_region(DT)
  
  # Sort by the sums, for easier comparison
  setkey(region_sums, "slg")
  
  expect_equal(region_sums[, slg], 
               c(1:8, c(10, 100, 1000, 10000)))
})



