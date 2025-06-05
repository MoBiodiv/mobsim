test_that("fails correctly with out-of-range values", {
  expect_error(
    sampling_grids(
      n_quadrats = 2L,
      xmin = 0.1,
      xmax = 1, # xmin is larger than x0
      ymin = 0,
      ymax = 1,
      x0 = 0,
      y0 = 0,
      delta_x = 0.1,
      delta_y = 0.1,
      quadrat_size = 0.01
    )
  )

  expect_error(
    sampling_grids(
      n_quadrats = 2L,
      xmin = 0,
      xmax = 1,
      ymin = 0.1,
      ymax = 1, # ymin is larger than y0
      x0 = 0,
      y0 = 0,
      delta_x = 0.1,
      delta_y = 0.1,
      quadrat_size = 0.01
    )
  )

  expect_warning(
    sampling_grids(
      n_quadrats = 10L,
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      x0 = 0.4,
      y0 = 0.4,
      delta_x = 0.01,
      delta_y = 0.01,
      quadrat_size = 0.1
    )
  )
})

test_that("classes are correct", {
  expect_type(
    sampling_grids(
      n_quadrats = 2L,
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      x0 = 0,
      y0 = 0,
      delta_x = 0.1,
      delta_y = 0.1,
      quadrat_size = 0.01
    ),
    "list"
  )
  expect_s3_class(
    sampling_grids(
      n_quadrats = 2L,
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      x0 = 0,
      y0 = 0,
      delta_x = 0.1,
      delta_y = 0.1,
      quadrat_size = 0.01
    ),
    "data.frame"
  )
})

test_that("dimensions are correct", {
  n_quadrats = 5L
  xy_dat <- sampling_grids(
    n_quadrats = n_quadrats,
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1,
    x0 = 0,
    y0 = 0,
    delta_x = 0.1,
    delta_y = 0.1,
    quadrat_size = 0.01
  )

  expect_equal(dim(xy_dat), c(n_quadrats, 2L))
})

test_that("samples are within range", {
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  quadrat_size = 0.01

  xy_dat <- sampling_grids(
    n_quadrats = 81L,
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    x0 = 0,
    y0 = 0,
    delta_x = 0.1,
    delta_y = 0.1,
    quadrat_size = quadrat_size
  )

  expect_gte(min(xy_dat$x), xmin)
  expect_lte(max(xy_dat$x) + quadrat_size, xmax)
  expect_gte(min(xy_dat$y), ymin)
  expect_lte(max(xy_dat$y) + quadrat_size, ymax)
})
