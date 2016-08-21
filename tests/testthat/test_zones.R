context("zone-creating functions")

test_that("dist_to_knn: returns correct order", {
  coords <- matrix(c(c(0, 0),
                     c(1, 0), 
                     c(4, 0),
                     c(1, 2),
                     c(-0.5, 2)),
                   ncol = 2, byrow = TRUE)
  m <- as.matrix(dist(coords, diag = T, upper = T))
  
  true_nns <- matrix(c(c(1, 2, 5, 4, 3),
                       c(2, 1, 4, 5, 3), 
                       c(3, 2, 4, 1, 5),
                       c(4, 5, 2, 1, 3),
                       c(5, 4, 1, 2, 3)),
                     ncol = 5, byrow = TRUE)
  
  expect_equal(unname(dist_to_knn(m)), true_nns)
})

test_that("k_nearest_neighbors: returns correct order", {
  coords <- matrix(c(c(0, 0),
                     c(1, 0), 
                     c(4, 0),
                     c(1, 2),
                     c(-0.5, 2)),
                   ncol = 2, byrow = TRUE)
  
  true_nns <- matrix(c(c(1, 2, 5, 4, 3),
                       c(2, 1, 4, 5, 3), 
                       c(3, 2, 4, 1, 5),
                       c(4, 5, 2, 1, 3),
                       c(5, 4, 1, 2, 3)),
                     ncol = 5, byrow = TRUE)
  
  expect_equal(unname(coords_to_knn(coords)), true_nns)
})

# test_that("k_nearest_neighbors: returns correct order", {
#   x <- matrix(c(c(-0.70,  1.01,  1.13, -0.11),
#                 c(-0.18,  0.82,  0.81, -0.76), 
#                 c(-0.14, -0.21,  0.33, -0.35),
#                 c(0.28,   0.65,  1.02,  0.35),
#                 c(0.40,   0.18, -0.59,  0.79)),
#               ncol = 4, byrow = TRUE)
#   
#   nn <- matrix(c(c(1, 2, 4, 3, 5),
#                  c(2, 1, 3, 4, 5), 
#                  c(3, 2, 4, 1, 5),
#                  c(4, 1, 2, 3, 5),
#                  c(5, 3, 4, 2, 1)),
#                ncol = 5, byrow = TRUE)
#   
#   expect_equal(unname(k_nearest_neighbors(x)), nn)
# })


test_that("closest_subsets: returns correct sets", {
  expres <- sets::set(sets::as.set(1L),
                      sets::as.set(1:2),
                      sets::as.set(1:3),
                      sets::as.set(1:4))
  expect_equal(closest_subsets(1:4), expres)
})

test_that("knn_zones: returns correct sets", {
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
  regs <- lapply(regs, function(x) unlist(as.list(x)))
  expect_equal(knn_zones(nn[, 1:2]), regs)
})

# Flexible zone shape (Tango 2005) -------------------------------------------

test_that("flexible_zones: works", {
  A <- matrix(c(0,1,0,0,0,0,
                1,0,1,0,0,0,
                0,1,0,0,0,0,
                0,0,0,0,1,0,
                0,0,0,1,0,0,
                0,0,0,0,0,0), 
              nrow = 6, byrow = TRUE) == 1
  kn <- matrix(as.integer(
               c(1,2,3,4,5,6,
                 2,1,3,4,5,6,
                 3,2,1,4,5,6,
                 4,5,1,6,3,2,
                 5,4,6,1,3,2,
                 6,5,4,1,3,2)),
               nrow = 6, byrow = TRUE)
  tru <- sets::set(sets::set(1L),
                         sets::set(2L),
                         sets::set(3L),
                         sets::set(4L),
                         sets::set(5L),
                         sets::set(6L),
                         sets::set(1L, 2L),
                         sets::set(2L, 3L),
                         sets::set(4L, 5L),
                         sets::set(1L, 2L, 3L))
  tru <- lapply(tru, FUN = function(x) unlist(as.list(x)))
  expect_equal(flexible_zones(kn, A), tru)
})

test_that("connected_neighbors: works", {
  A <- matrix(c(0,1,0,0,0,0,
                1,0,1,0,0,0,
                0,1,0,0,0,0,
                0,0,0,0,1,0,
                0,0,0,1,0,0,
                0,0,0,0,0,0), 
              nrow = 6, byrow = TRUE)
  A <- A == 1
  
  expect_equal(connected_neighbors(1:6, A), 
               sets::set(sets::set(1L), 
                         sets::set(1L, 2L),
                         sets::set(1L, 2L, 3L)))
  expect_equal(connected_neighbors(c(2:6, 1L), A), 
               sets::set(sets::set(2L), 
                         sets::set(1L, 2L),
                         sets::set(2L, 3L),
                         sets::set(1L, 2L, 3L)))
  expect_equal(connected_neighbors(c(3:6, 1:2), A), 
               sets::set(sets::set(3L), 
                         sets::set(2L, 3L),
                         sets::set(1L, 2L, 3L)))
  expect_equal(connected_neighbors(c(4:6, 1:3), A), 
               sets::set(sets::set(4L), 
                         sets::set(4L, 5L)))
  expect_equal(connected_neighbors(c(5:6, 1:4), A), 
               sets::set(sets::set(5L), 
                         sets::set(4L, 5L)))
  expect_equal(connected_neighbors(c(6L, 1:5), A), 
               sets::set(sets::set(6L)))
})

test_that("if_connected: works", {
  A <- matrix(c(0,1,0,0,0,
                1,0,1,0,0,
                0,1,0,0,0,
                0,0,0,0,1,
                0,0,0,1,0), 
              nrow = 5, byrow = TRUE)
  A <- A == 1
  
  expect_equal(if_connected(sets::set(2L), 1L, A), sets::set(1L, 2L))
  expect_equal(if_connected(sets::set(2L, 3L), 1L, A), sets::set(1L, 2L, 3L))
  expect_equal(if_connected(sets::set(4L), 1L, A), sets::set())
  expect_equal(if_connected(sets::set(2L, 4L), 1L, A), sets::set())
})


test_that("is_connected: works", {
  A <- matrix(c(0,1,0,0,0,
                1,0,1,0,0,
                0,1,0,0,0,
                0,0,0,0,1,
                0,0,0,1,0), 
              nrow = 5, byrow = TRUE)
  A <- A == 1
  expect_true(is_connected(sets::set(2L), 1L, A))
  expect_true(is_connected(sets::set(2L, 3L), 1L, A))
  expect_false(is_connected(sets::set(4L), 1L, A))
  expect_false(is_connected(sets::set(2L, 4L), 1L, A))
})



test_that("connected_to: works", {
  A <- matrix(c(0,1,0,0,0,
                1,0,1,0,0,
                0,1,0,0,0,
                0,0,0,0,1,
                0,0,0,1,0), 
              nrow = 5, byrow = TRUE)
  A <- A == 1
  z0a <- sets::as.set(1L)
  z1a <- sets::as.set(2L)
  actual_a <- connected_to(z0a, z1a, A)
  
  z0b <- sets::as.set(1L)
  z1b <- sets::set(4L, 5L)
  actual_b <- connected_to(z0b, z1b, A)
  
  z0c <- sets::as.set(2L)
  z1c <- sets::set(1L, 3L)
  actual_c <- connected_to(z0c, z1c, A)
  
  expect_equal(actual_a, sets::set(2L))
  expect_equal(actual_b, sets::set())
  expect_equal(actual_c, sets::set(1L, 3L))
})