test_that("snapshot test with clustering", {
  mpcoords <- create_motherpoints(s_local = 100L, seed = 42L)
  expect_snapshot(drift_x_motherpoints(mpcoords, drift = 0.1))
})

test_that("snapshot test without clustering", {
  mpcoords <- create_motherpoints(
    s_local = 100L,
    n_motherpoints_range = 0:1,
    seed = 42L
  )
  expect_snapshot(drift_x_motherpoints(mpcoords, drift = 0.1))
})

test_that("snapshot test with clustering", {
  mpcoords <- create_motherpoints(s_local = 100L, seed = 42L)
  expect_snapshot(drift_y_motherpoints(mpcoords, drift = 0.1))
})

test_that("snapshot test without clustering", {
  mpcoords <- create_motherpoints(
    s_local = 100L,
    n_motherpoints_range = 0:1,
    seed = 42L
  )
  expect_snapshot(drift_y_motherpoints(mpcoords, drift = 0.1))
})
