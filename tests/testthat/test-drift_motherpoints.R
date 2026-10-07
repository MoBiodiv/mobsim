## invariants ----
test_that("drift_x_motherpoints() shifts all x by exactly drift", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  d <- 0.15
  result <- drift_x_motherpoints(mpcoords, drift = d)
  expect_equal(
    mapply(
      FUN = function(r, o) all(r == o + d),
      result$xmother,
      mpcoords$xmother
    ),
    rep(TRUE, length(mpcoords$xmother))
  )
})

test_that("drift_y_motherpoints() shifts all y by exactly drift", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  d <- -0.3
  result <- drift_y_motherpoints(mpcoords, drift = d)
  expect_equal(
    mapply(
      FUN = function(r, o) all(r == o + d),
      result$ymother,
      mpcoords$ymother
    ),
    rep(TRUE, length(mpcoords$ymother))
  )
})

test_that("drift_x_motherpoints() does not change y coordinates", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  result <- drift_x_motherpoints(mpcoords, drift = 0.1)
  expect_equal(result$ymother, mpcoords$ymother)
})

test_that("drift_y_motherpoints() does not change x coordinates", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  result <- drift_y_motherpoints(mpcoords, drift = 0.1)
  expect_equal(result$xmother, mpcoords$xmother)
})

test_that("drift_x_motherpoints() preserves mother_map class", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  result <- drift_x_motherpoints(mpcoords, drift = 0.1)
  expect_s3_class(result, "mother_map")
})

test_that("drift_y_motherpoints() preserves mother_map class", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  result <- drift_y_motherpoints(mpcoords, drift = 0.1)
  expect_s3_class(result, "mother_map")
})

## edge cases ----
test_that("drift_x_motherpoints() errors when drift is not provided", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  expect_error(drift_x_motherpoints(mpcoords))
})

test_that("drift_y_motherpoints() errors when drift is not provided", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  expect_error(drift_y_motherpoints(mpcoords))
})

test_that("drift_x_motherpoints() with zero drift leaves x unchanged", {
  mpcoords <- create_motherpoints(s_local = 10L, seed = 42L)
  result <- drift_x_motherpoints(mpcoords, drift = 0)
  expect_equal(result$xmother, mpcoords$xmother)
})

## snapshots ----
test_that("snapshot test drift_x with clustering", {
  mpcoords <- create_motherpoints(s_local = 100L, seed = 42L)
  expect_snapshot(drift_x_motherpoints(mpcoords, drift = 0.1))
})

test_that("snapshot test drift_x without clustering", {
  mpcoords <- create_motherpoints(
    s_local = 100L,
    n_motherpoints_range = 0:1,
    seed = 42L
  )
  expect_snapshot(drift_x_motherpoints(mpcoords, drift = 0.1))
})

test_that("snapshot test drift_y with clustering", {
  mpcoords <- create_motherpoints(s_local = 100L, seed = 42L)
  expect_snapshot(drift_y_motherpoints(mpcoords, drift = 0.1))
})

test_that("snapshot test drift_y without clustering", {
  mpcoords <- create_motherpoints(
    s_local = 100L,
    n_motherpoints_range = 0:1,
    seed = 42L
  )
  expect_snapshot(drift_y_motherpoints(mpcoords, drift = 0.1))
})
