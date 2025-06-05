test_that("fails correctly", {
  sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)

  expect_error(sample_quadrats(data.frame(c(1:5), c(1:5))))
  expect_error(sample_quadrats(sim_com1, n_quadrats = 10.5))
  expect_error(sample_quadrats(sim_com1, quadrat_area = 0.01, method = "foo"))
})

test_that("classes are as expected", {
  comm_mat1 <- sample_quadrats(
    comm = sim_poisson_community(s_pool = 5L, n_sim = 50L),
    plot = FALSE,
    n_quadrats = 2L,
    method = "random",
    avoid_overlap = FALSE
  )

  expect_type(comm_mat1, "list")
  expect_s3_class(comm_mat1[[1]], "data.frame")
  expect_s3_class(comm_mat1[[2]], "data.frame")
})

test_that("dimensions are as expected", {
  S = 5L
  sim_com1 <- sim_poisson_community(
    s_pool = S,
    n_sim = 50L,
    fix_s_sim = TRUE
  )
  n_quadrats <- 10L
  comm_mat <- sample_quadrats(
    sim_com1,
    n_quadrats = n_quadrats,
    method = "random",
    plot = FALSE,
    avoid_overlap = TRUE
  )

  expect_vector(comm_mat, size = 2L)
  expect_equal(dim(comm_mat[[1]]), c(n_quadrats, S))
  expect_equal(dim(comm_mat[[2]]), c(n_quadrats, 2L))
})

test_that("default argument values are as expected", {
  comm <- sim_thomas_community(
    s_pool = 10L,
    n_sim = 100L,
    sad_type = "lnorm",
    sad_coef = list(cv_abund = 1),
    fix_s_sim = TRUE,
    sigma = 1,
    mother_points = 0L,
    xrange = c(0, 1),
    yrange = c(0, 1),
    seed = 42L
  )
  expect_snapshot_output(sample_quadrats(comm = comm, seed = 42L))
})

test_that("calls the correct subfunction if method is random, with overlap", {
  sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
  quadrat_area <- 0.01
  expect_equal(
    as.matrix(
      suppressWarnings(
        sample_quadrats(
          comm = sim_com1,
          seed = 42L,
          plot = FALSE,
          n_quadrats = 10L,
          quadrat_area = quadrat_area,
          method = "random",
          avoid_overlap = FALSE
        )$xy_dat
      ),
      rownames.force = FALSE
    ),
    as.matrix(
      suppressWarnings(
        sampling_random_overlap(
          xmin = sim_com1$x_min_max[1L],
          xmax = sim_com1$x_min_max[2L] - sqrt(quadrat_area),
          ymin = sim_com1$y_min_max[1L],
          ymax = sim_com1$y_min_max[2L] - sqrt(quadrat_area),
          min_dist = sqrt(2 * quadrat_area),
          n_quadrats = 10L,
          seed = 42L
        )
      )
    )
  )
})

test_that("calls the correct subfunction if method is random, without overlap, with spatstat", {
  testthat::skip_if_not(requireNamespace("spatstat.random", quietly = TRUE))
  sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
  quadrat_area <- 0.01
  expect_equal(
    as.matrix(
      sample_quadrats(
        comm = sim_com1,
        seed = 42L,
        plot = FALSE,
        n_quadrats = 10L,
        quadrat_area = quadrat_area,
        method = "random",
        avoid_overlap = TRUE
      )$xy_dat,
      rownames.force = FALSE
    ),
    as.matrix(
      sampling_random_spatstat(
        xmin = sim_com1$x_min_max[1L],
        xmax = sim_com1$x_min_max[2L] - sqrt(quadrat_area),
        ymin = sim_com1$y_min_max[1L],
        ymax = sim_com1$y_min_max[2L] - sqrt(quadrat_area),
        min_dist = sqrt(2 * quadrat_area),
        n_quadrats = 10L,
        seed = 42L
      )
    )
  )
})

test_that("calls the correct subfunction if method is random, without overlap, without spatstat", {
  mockery::stub(
    where = sample_quadrats,
    what = "requireNamespace",
    how = FALSE,
    depth = 1L
  )

  sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
  quadrat_area <- 0.01
  expect_equal(
    as.matrix(
      sample_quadrats(
        comm = sim_com1,
        seed = 42L,
        plot = FALSE,
        n_quadrats = 10L,
        quadrat_area = quadrat_area,
        method = "random",
        avoid_overlap = TRUE
      )$xy_dat,
      rownames.force = FALSE
    ),
    as.matrix(
      sampling_random_bruteforce(
        xmin = sim_com1$x_min_max[1L],
        xmax = sim_com1$x_min_max[2L] - sqrt(quadrat_area),
        ymin = sim_com1$y_min_max[1L],
        ymax = sim_com1$y_min_max[2L] - sqrt(quadrat_area),
        min_dist = sqrt(2 * quadrat_area),
        n_quadrats = 10L,
        seed = 42L
      )
    )
  )
})

test_that("calls the correct subfunction if method is transects", {
  sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
  quadrat_area <- 0.01
  expect_equal(
    as.matrix(
      sample_quadrats(
        comm = sim_com1,
        plot = FALSE,
        n_quadrats = 5L,
        quadrat_area = quadrat_area,
        method = "transect"
      )$xy_dat,
      rownames.force = FALSE
    ),
    as.matrix(
      sampling_transects(
        xmin = sim_com1$x_min_max[1L],
        xmax = sim_com1$x_min_max[2L] - sqrt(quadrat_area),
        ymin = sim_com1$y_min_max[1L],
        ymax = sim_com1$y_min_max[2L] - sqrt(quadrat_area),
        x0 = 0,
        y0 = 0,
        delta_x = 0.1,
        delta_y = 0.1,
        quadrat_size = sqrt(quadrat_area),
        n_quadrats = 5L
      )
    )
  )
})
test_that("calls the correct subfunction if method is grid", {
  sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
  quadrat_area <- 0.01
  expect_equal(
    as.matrix(
      sample_quadrats(
        sim_com1,
        plot = FALSE,
        n_quadrats = 5L,
        quadrat_area = quadrat_area,
        method = "grid"
      )$xy_dat,
      rownames.force = FALSE
    ),
    as.matrix(
      sampling_grids(
        xmin = sim_com1$x_min_max[1L],
        xmax = sim_com1$x_min_max[2L] - sqrt(quadrat_area),
        ymin = sim_com1$y_min_max[1L],
        ymax = sim_com1$y_min_max[2L] - sqrt(quadrat_area),
        x0 = 0,
        y0 = 0,
        delta_x = 0.1,
        delta_y = 0.1,
        quadrat_size = sqrt(quadrat_area),
        n_quadrats = 5L
      ),
      rownames.force = FALSE
    )
  )
})

test_that("calls the correct subfunction if edge case where n_quadrats == 1", {
  sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
  quadrat_area <- 0.01
  expect_equal(
    as.numeric(
      sample_quadrats(
        comm = sim_com1,
        seed = 42L,
        plot = FALSE,
        n_quadrats = 1L,
        quadrat_area = quadrat_area,
        method = "random",
        avoid_overlap = TRUE
      )$xy_dat
    ),
    as.numeric(
      sampling_one_quadrat(
        xmin = sim_com1$x_min_max[1L],
        xmax = sim_com1$x_min_max[2L] - sqrt(quadrat_area),
        ymin = sim_com1$y_min_max[1L],
        ymax = sim_com1$y_min_max[2L] - sqrt(quadrat_area),
        seed = 42L
      )
    )
  )
})

test_that("plotting works", {
  sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
  expect_snapshot_output(
    sample_quadrats(
      comm = sim_com1,
      n_quadrats = 15L,
      quadrat_area = 0.001,
      seed = 42L,
      plot = TRUE,
      method = "random",
      avoid_overlap = TRUE
    )
  )
  expect_snapshot_output(
    sample_quadrats(
      comm = sim_com1,
      n_quadrats = 25L,
      quadrat_area = 0.001,
      seed = 42L,
      plot = TRUE,
      method = "grid"
    )
  )
})
