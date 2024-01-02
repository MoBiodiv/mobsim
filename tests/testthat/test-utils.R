# tests utility functions
## torusify() ----
test_that("torusify() - mother_map method - default snapshot 0-1 range", {
   mpcoords <- create_motherpoints(s_local = 20L, seed = 42L)
   mpcoordsJ <- jitter_motherpoints(mpcoords, sd = 3, seed = 42L)
   expect_snapshot(torusify(mpcoordsJ))
})

test_that("torusify() - mother_map method - wild range snapshot", {
   mpcoords <- create_motherpoints(s_local = 20L, xrange = c(8, 10), yrange = c(-1, 1), seed = 42L)
   mpcoordsJ <- jitter_motherpoints(mpcoords, sd = 3, seed = 42L)
   expect_snapshot(torusify(mpcoordsJ))
})

test_that("torusify() - community method - wild range snapshot", {
   simdat <- mobsim::sim_thomas_community(10L, 10L, xrange = c(8, 10), yrange = c(-1, 1), seed = 42L)
   simdatJ <- jitter_species(simdat, sd = 3, seed = 42L)
   expect_snapshot(torusify(simdatJ))
})

test_that("torusify() - data.frame method - wild range snapshot", {
   set.seed(42L)
   x_min_max <- c(8, 10)
   y_min_max <- c(-1, 1)
   mpcoords <- data.frame(
      xmother = stats::runif(20L, min = x_min_max[1L] - 1, max = x_min_max[2L] + 1),
      ymother = stats::runif(20L, min = y_min_max[1L] - 1, max = y_min_max[2L] + 1)
   )
   expect_snapshot(torusify(mpcoords, x_min_max, y_min_max))

})

test_that("torudify() - numeric method - wild range snapshot", {
   set.seed(42L)
   mpcoords <- stats::runif(10L, min = -2, max = 2)
   expect_snapshot(torusify(mpcoords, range = c(0, 1)))
})


## community2ppp() ----
test_that("community2ppp() default behaviour snapshot", {
   simdat <- mobsim::sim_thomas_coords(11L:100L)
   expect_snapshot(community2ppp(simdat))
})

## create_random_ID() ----
test_that("create_random_ID() default behaviour", {
   expect_equal(create_random_ID(seed = 42L), "QEAYJ1252R")
})

test_that("create_random_ID() is not affected by global seed", {
   set.seed(84L)
   expect_false(create_random_ID() == create_random_ID(seed = 84L))
})

test_that("create_random_ID() does not affect global seed", {
   set.seed(84L)
   reference84 <- .Random.seed
   create_random_ID(seed = 42L)
   expect_equal(reference84, .Random.seed)
})
