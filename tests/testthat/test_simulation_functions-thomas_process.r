# sim_thomas_coords, sim_thomas_community, rThomas_rcpp,
test_that("sim_thomas_community() - default argument values are as expected", {
  expect_snapshot(sim_thomas_community(s_pool = 3L, n_sim = 100L, seed = 42L))
})

test_that("sim_thomas_coords() - default argument values are as expected", {
  expect_snapshot(sim_thomas_coords(abund_vec = 1:4, seed = 42L))
})

test_that("rThomas_rcpp() - default argument values are as expected", {
  set.seed(42L)
  expect_snapshot(rThomas_rcpp(
    n_points = 100L,
    n_mother_points = 2L,
    xmother = NA,
    ymother = NA,
    sigma = 0.2
  ))
  expect_snapshot(rThomas_rcpp(
    n_points = 100L,
    n_mother_points = 2L,
    xmother = 0.5,
    ymother = 0.5,
    sigma = 0.2
  ))
  expect_snapshot({
    rThomas_rcpp(
      n_points = 100L,
      n_mother_points = 2L,
      xmother = c(0.2, 0.8),
      ymother = c(0.2, 0.8),
      sigma = 0.2
    )
  })
})

test_that("sim_thomas_community() calls sim_sad() correctly", {
  s_pool <- 10L
  n_sim <- 100L

  sad1 <- sim_sad(s_pool = s_pool, n_sim = n_sim, fix_s_sim = TRUE, seed = 42L)
  sim1 <- sim_thomas_community(
    s_pool = s_pool,
    n_sim = n_sim,
    fix_s_sim = TRUE,
    seed = 42L
  )

  expect_equal(sum(sad1), nrow(sim1$census))
  expect_equal(length(sad1), length(table(sim1$census$species)))
})

test_that("sim_thomas_coords() - handles wrong mother_points parametres", {
  expect_error({
    # mother_points is negative
    sim_thomas_coords(abund_vec = sad1, mother_points = -1)
  })
  expect_error({
    # mother_points contains negative values
    sim_thomas_coords(
      abund_vec = sad1,
      mother_points = sample(size = 100, x = -2:2, replace = TRUE)
    )
  })
  expect_error({
    # mother_points is too short
    sim_thomas_coords(abund_vec = sad1, mother_points = rep(1, 4))
  })
  expect_error({
    # mother_points is too long
    sim_thomas_coords(abund_vec = sad1, mother_points = rep(1, 1000))
  })
})

test_that("sim_thomas_coords() - handles wrong xmother and ymother parametres", {
  sad2 <- sim_sad(s_pool = 2, n_sim = 100)
  xmother <- lapply(1:length(sad2), function(x) runif(2, 0, 1))

  expect_error(sim_thomas_coords(
    sad2,
    xmother = xmother,
    ymother = list(c(0.2, NA), c(0.3, NA))
  ))
  expect_error(sim_thomas_coords(
    sad2,
    xmother = xmother,
    ymother = list(NA, c(0.2, 0.3))
  ))
  expect_error(sim_thomas_coords(
    sad2,
    xmother = xmother,
    ymother = list(2, c(0.2, 0.3))
  ))
})

test_that("sim_thomas_coords() - handles wrong xrange and yrange parametres", {
  expect_error(sim_thomas_coords(
    sad2,
    xmother = xmother,
    ymother = xmother,
    xrange = c(0, 0.1)
  )) # xmother and ymother outside of range
  expect_error(sim_thomas_coords(
    sad2,
    xmother = xmother,
    ymother = xmother,
    xrange = data.frame(c(0, 1), c(0, 0.1))
  )) # xrange and y range have different class
  expect_error(sim_thomas_coords(
    sad2,
    xmother = xmother,
    ymother = xmother,
    xrange = data.frame(c(0, 1, 2), c(0, 1, 2)),
    yrange = data.frame(c(0, 1, 2), c(0, 1, 2))
  )) # xrange and yrange have too many columns
  expect_error(sim_thomas_coords(
    sad2,
    xmother = xmother,
    ymother = xmother,
    xrange = data.frame(c(0, 1), c(0, 1), c(0, 1)),
    yrange = data.frame(c(0, 1), c(0, 1), c(0, 1))
  )) # xrange and yrange have too many rows
})

test_that("sim_thomas_coords() - throws warnings when expected", {
  sad1 <- sim_sad(s_pool = 3L, n_sim = 10L)
  xmother <- lapply(1L:length(sad1), function(x) runif(1, 0, 1))
  ymother <- lapply(1L:length(sad1), function(x) runif(1, 0, 1))
  expect_warning(sim_thomas_coords(
    sad1,
    mother_points = 1,
    xmother = xmother,
    ymother = ymother
  ))
  expect_warning(sim_thomas_coords(sad1, mother_points = 1, cluster_points = 2))
})

test_that("rThomas_rcpp() - behaves as expected", {
  expect_s3_class(
    rThomas_rcpp(
      20,
      n_mother_points = 1,
      xmother = .5,
      ymother = .5,
      sigma = 0.2
    ),
    "data.frame"
  )
  expect_equal(
    nrow(rThomas_rcpp(
      20,
      n_mother_points = 1,
      xmother = .5,
      ymother = .5,
      sigma = 0.2
    )),
    20
  )
  expect_error(rThomas_rcpp(20, n_mother_points = 2, sigma = 0.2))
})
