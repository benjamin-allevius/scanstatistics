context("Region-creating functions")

test_that("k_nearest_neighbors: returns correct order", {
  x <- matrix(c(c(-0.70,  1.01,  1.13, -0.11),
                c(-0.18,  0.82,  0.81, -0.76), 
                c(-0.14, -0.21,  0.33, -0.35),
                c(0.28,   0.65,  1.02,  0.35),
                c(0.40,   0.18, -0.59,  0.79)),
              ncol = 4, byrow = TRUE)
  
  nn <- matrix(c(c(1, 2, 4, 3, 5),
                 c(2, 1, 3, 4, 5), 
                 c(3, 2, 4, 1, 5),
                 c(4, 1, 2, 3, 5),
                 c(5, 3, 4, 2, 1)),
               ncol = 5, byrow = TRUE)
  
  expect_equal(unname(k_nearest_neighbors(x)), nn)
})


test_that("closest_subsets: returns correct sets", {
  expres <- sets::set(sets::as.set(1L),
                      sets::as.set(1:2),
                      sets::as.set(1:3),
                      sets::as.set(1:4))
  expect_equal(closest_subsets(1:4), expres)
})

test_that("regions_upto_k: returns correct sets", {
  nn <- matrix(c(c(1L, 2L, 4L, 3L, 5L),
                 c(2L, 1L, 3L, 4L, 5L), 
                 c(3L, 2L, 4L, 1L, 5L),
                 c(4L, 1L, 2L, 3L, 5L),
                 c(5L, 3L, 4L, 2L, 1L)),
               ncol = 5, byrow = TRUE)
  regs <- sets::set(sets::as.set(1L),
                    sets::as.set(2L),
                    sets::as.set(3L),
                    sets::as.set(4L),
                    sets::as.set(5L),
                    sets::as.set(c(1L, 2L)),
                    sets::as.set(c(3L, 2L)),
                    sets::as.set(c(4L, 1L)),
                    sets::as.set(c(5L, 3L)))
  expect_equal(regions_upto_k(nn[, 1:2]), regs)
})


test_that("is_connected: works", {
  A <- matrix(c(0,1,0,0,0,
                1,0,1,0,0,
                0,1,0,0,0,
                0,0,0,0,1,
                0,0,0,1,0), 
              nrow = 5, byrow = TRUE)
  A <- A == 1
  connected_to <- pryr::partial(connected_to_full,
                                adjacency_matrix = A)
  expect_true(is_connected(sets::set(2L), 1L))
  expect_true(is_connected(sets::set(2L, 3L), 1L))
  expect_false(is_connected(sets::set(4L), 1L))
  expect_false(is_connected(sets::set(2L, 4L), 1L))
})



test_that("connected_to_full: works", {
  A <- matrix(c(0,1,0,0,0,
                1,0,1,0,0,
                0,1,0,0,0,
                0,0,0,0,1,
                0,0,0,1,0), 
              nrow = 5, byrow = TRUE)
  A <- A == 1
  z0a <- sets::as.set(1L)
  z1a <- sets::as.set(2L)
  actual_a <- connected_to_full(z0a, z1a, A)
  
  z0b <- sets::as.set(1L)
  z1b <- sets::set(4L, 5L)
  actual_b <- connected_to_full(z0b, z1b, A)
  
  z0c <- sets::as.set(2L)
  z1c <- sets::set(1L, 3L)
  actual_c <- connected_to_full(z0c, z1c, A)
  
  expect_equal(actual_a, sets::set(2L))
  expect_equal(actual_b, sets::set())
  expect_equal(actual_c, sets::set(1L, 3L))
})