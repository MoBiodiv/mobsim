## structure ----
test_that("create_motherpoints() output is class mother_map", {
  result <- create_motherpoints(s_local = 10L, seed = 1L)
  expect_s3_class(result, "mother_map")
  expect_true(all(is.element(
    c("n_mother_points", "xmother", "ymother", "x_min_max", "y_min_max"),
    names(result)
  )))
})

test_that("create_motherpoints() xmother and ymother are always lists", {
  # Default range — sapply would have given a list anyway, confirm list
  r1 <- create_motherpoints(s_local = 10L, seed = 1L)
  expect_type(r1$xmother, "list")
  expect_type(r1$ymother, "list")

  # Uniform count (n_motherpoints_range = 1L) — previously collapsed to vector
  r2 <- create_motherpoints(s_local = 10L, n_motherpoints_range = 1L, seed = 1L)
  expect_type(r2$xmother, "list")
  expect_type(r2$ymother, "list")

  # Uniform count > 1 — previously collapsed to matrix
  r3 <- create_motherpoints(s_local = 10L, n_motherpoints_range = 2L, seed = 1L)
  expect_type(r3$xmother, "list")
  expect_type(r3$ymother, "list")
})

test_that("create_motherpoints() list lengths match s_local", {
  s <- 15L
  result <- create_motherpoints(s_local = s, seed = 1L)
  expect_length(result$xmother, s)
  expect_length(result$ymother, s)
  expect_length(result$n_mother_points, s)
})

## invariants ----
test_that("create_motherpoints() coordinates stay within xrange/yrange", {
  xr <- c(2, 5)
  yr <- c(-1, 3)
  result <- create_motherpoints(
    s_local = 20L,
    xrange = xr,
    yrange = yr,
    seed = 1L
  )
  all_x <- unlist(result$xmother)
  all_y <- unlist(result$ymother)
  expect_true(all(all_x >= xr[[1L]] & all_x <= xr[[2L]]))
  expect_true(all(all_y >= yr[[1L]] & all_y <= yr[[2L]]))
})

test_that("create_motherpoints() n_mother_points values are within n_motherpoints_range", {
  rng <- 1L:5L
  result <- create_motherpoints(
    s_local = 50L,
    n_motherpoints_range = rng,
    seed = 1L
  )
  expect_true(all(is.element(result$n_mother_points, rng)))
})

## reproducibility ----
test_that("create_motherpoints() seed produces identical results", {
  r1 <- create_motherpoints(s_local = 20L, seed = 99L)
  r2 <- create_motherpoints(s_local = 20L, seed = 99L)
  expect_equal(r1$n_mother_points, r2$n_mother_points)
  expect_equal(r1$xmother, r2$xmother)
  expect_equal(r1$ymother, r2$ymother)
})

test_that("create_motherpoints() stores xrange and yrange correctly", {
  xr <- c(2, 5)
  yr <- c(-1, 3)
  result <- create_motherpoints(
    s_local = 5L,
    xrange = xr,
    yrange = yr,
    seed = 1L
  )
  expect_equal(result$x_min_max, xr)
  expect_equal(result$y_min_max, yr)
})

## snapshots ----
test_that("snapshot test with clustering", {
  expect_snapshot(create_motherpoints(s_local = 100L, seed = 42L))
})

test_that("snapshot test with clustering - no range", {
  expect_snapshot(create_motherpoints(
    s_local = 100L,
    n_motherpoints_range = 1L,
    seed = 42L
  ))
})

test_that("snapshot test without clustering", {
  expect_snapshot(create_motherpoints(
    s_local = 100L,
    n_motherpoints_range = 0L:1L,
    seed = 42L
  ))
})

test_that("snapshot test plotting method for mother_map class", {
  mpcoords <- create_motherpoints(
    s_local = 100L,
    n_motherpoints_range = 0L:1L,
    seed = 42L
  )
  expect_snapshot_output(plot(mpcoords))
})
