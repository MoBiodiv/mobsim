## structure ----
test_that("jitter_motherpoints() returns a mother_map object", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  result <- jitter_motherpoints(mpcoords, seed = 1L)
  expect_s3_class(result, "mother_map")
  expect_true(all(is.element(
    c("n_mother_points", "xmother", "ymother", "x_min_max", "y_min_max"),
    names(result)
  )))
})

test_that("jitter_motherpoints() preserves list lengths", {
  s <- 15L
  mpcoords <- create_motherpoints(s_local = s, seed = 42L)
  result <- jitter_motherpoints(mpcoords, seed = 1L)
  expect_length(result$xmother, s)
  expect_length(result$ymother, s)
  expect_length(result$n_mother_points, s)
})

test_that("jitter_motherpoints() preserves n_mother_points", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  result <- jitter_motherpoints(mpcoords, seed = 1L)
  expect_equal(result$n_mother_points, mpcoords$n_mother_points)
})

## invariants ----
test_that("jitter_motherpoints() coordinates actually change", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  result <- jitter_motherpoints(mpcoords, sd = 0.1, seed = 1L)
  expect_false(isTRUE(all.equal(
    unlist(result$xmother),
    unlist(mpcoords$xmother)
  )))
})

test_that("jitter_motherpoints() species with 0 mother points are unaffected", {
  mpcoords <- create_motherpoints(
    s_local = 20L,
    n_motherpoints_range = 0L:3L,
    seed = 42L
  )
  zero_idx <- which(mpcoords$n_mother_points == 0L)
  if (length(zero_idx) > 0L) {
    result <- jitter_motherpoints(mpcoords, sd = 1, seed = 1L)
    for (i in zero_idx) {
      expect_length(result$xmother[[i]], 0L)
      expect_length(result$ymother[[i]], 0L)
    }
  }
})

## reproducibility ----
test_that("jitter_motherpoints() seed produces identical results", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  r1 <- jitter_motherpoints(mpcoords, seed = 7L)
  r2 <- jitter_motherpoints(mpcoords, seed = 7L)
  expect_equal(r1$xmother, r2$xmother)
  expect_equal(r1$ymother, r2$ymother)
})

## snapshots ----
test_that("snapshot test with clustering", {
  mpcoords <- create_motherpoints(s_local = 100L, seed = 42L)
  expect_snapshot(jitter_motherpoints(mpcoords, seed = 42L))
})

test_that("snapshot test without clustering", {
  mpcoords <- create_motherpoints(
    s_local = 100L,
    n_motherpoints_range = 0:1,
    seed = 42L
  )
  expect_snapshot(jitter_motherpoints(mpcoords, seed = 42L))
})
