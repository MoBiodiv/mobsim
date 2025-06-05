test_that("sim_sad() - correct assertions", {
  expect_error(sim_sad(n_sim = 100L, sad_type = "lnorm"))
  expect_error(
    sim_sad(s_pool = NA, n_sim = 100L, sad_type = "lnorm"),
    "The argument s_pool is mandatory for the selected sad and has to be a positive integer number."
  )
  expect_error(
    sim_sad(s_pool = 3.5, n_sim = 100L, sad_type = "lnorm"),
    "s_pool has to be a positive integer number"
  )
  expect_error(
    sim_sad(s_pool = 0L, n_sim = 100L, sad_type = "lnorm"),
    "s_pool has to be a positive integer number"
  )

  expect_error(sim_sad(s_pool = 3L, sad_type = "lnorm"))
  expect_error(sim_sad(s_pool = 3L, n_sim = NA, sad_type = "lnorm"))
  expect_error(
    sim_sad(s_pool = 3L, n_sim = 100.5, sad_type = "lnorm"),
    "n_sim has to be a positive integer number"
  )
  expect_error(
    sim_sad(s_pool = 3L, n_sim = 0L, sad_type = "lnorm"),
    "n_sim has to be a positive integer number"
  )

  expect_error(
    sim_sad(s_pool = 3L, n_sim = 100L, sad_type = "lnorm", sad_coef = 0.1),
    "coef must be a named list!"
  )
  expect_error(sim_sad(s_pool = 3L, n_sim = 100L, sad_type = "wrong_sad_type"))

  expect_error(sim_sad(sad_type = "bs", sad_coef = list(S = 4L)))

  expect_warning(sim_sad(
    s_pool = 3L,
    sad_type = "bs",
    sad_coef = list(S = 4L, N = 10L)
  ))
  expect_warning(sim_sad(
    n_sim = 10L,
    sad_type = "bs",
    sad_coef = list(S = 4L, N = 10L)
  ))
})

test_that("default argument values are as expected", {
  expect_snapshot(sim_sad(s_pool = 3L, n_sim = 100L, seed = 42L))
})

test_that("sim_sad() - results are as expected", {
  expect_equal(names(sim_sad(s_pool = 3L, n_sim = 10L))[1L], "species_1")
  n_sim <- 10L
  s_pool <- 3L

  sadNorm <- sim_sad(
    s_pool = s_pool,
    n_sim = n_sim,
    sad_type = "lnorm",
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadNorm), n_sim)
  expect_equal(length(sadNorm), s_pool)

  sadBs <- sim_sad(
    sad_type = "bs",
    sad_coef = list(S = s_pool, N = n_sim),
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadBs), n_sim)
  expect_equal(length(sadBs), s_pool)

  sadGamma <- sim_sad(
    s_pool = s_pool,
    n_sim = n_sim,
    sad_type = "gamma",
    sad_coef = list(shape = 1, scale = 1),
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadGamma), n_sim)
  expect_equal(length(sadGamma), s_pool)

  sadGeom <- sim_sad(
    s_pool = s_pool,
    n_sim = n_sim,
    sad_type = "geom",
    sad_coef = list(prob = 0.5),
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadGeom), n_sim)
  expect_equal(length(sadGeom), s_pool)

  sadLs <- sim_sad(
    sad_type = "ls",
    sad_coef = list(N = n_sim, alpha = 3),
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadLs), n_sim)

  # expect_equal(sum(sim_sad(sad_type = "mzsm", sad_coef = list(n = 3L, J = n_sim, theta = 4))), n_sim)

  sadNbinom <- sim_sad(
    s_pool = s_pool,
    n_sim = n_sim,
    sad_type = "nbinom",
    sad_coef = list(size = 100L, prob = 0.5),
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadNbinom), n_sim)
  expect_equal(length(sadNbinom), s_pool)

  sadPareto <- sim_sad(
    s_pool = s_pool,
    n_sim = n_sim,
    sad_type = "pareto",
    sad_coef = list(shape = 2),
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadPareto), n_sim)
  expect_equal(length(sadPareto), s_pool)

  sadPoilog <- sim_sad(
    s_pool = s_pool,
    n_sim = n_sim,
    sad_type = "poilog",
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadPoilog), n_sim)
  expect_equal(length(sadPoilog), s_pool)

  sadPower <- sim_sad(
    s_pool = s_pool,
    n_sim = n_sim,
    sad_type = "power",
    sad_coef = list(s = 2),
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadPower), n_sim)
  expect_equal(length(sadPower), s_pool)

  # expect_equal(sum(sim_sad(s_pool = 3L, n_sim = n_sim, sad_type = "powbend", sad_coef = list(s = 2, omega = 3))), n_sim)

  sadWeibull <- sim_sad(
    s_pool = s_pool,
    n_sim = n_sim,
    sad_type = "weibull",
    sad_coef = list(shape = 2),
    fix_s_sim = TRUE
  )
  expect_equal(sum(sadWeibull), n_sim)
  expect_equal(length(sadWeibull), s_pool)
})

test_that("sim_sad() - edge case s_pool = 1", {
  n_sim <- 100L
  s_pool <- 1L
  sad1lnorm <- sim_sad(
    s_pool = s_pool,
    n_sim = n_sim,
    sad_type = "lnorm",
    sad_coef = list(cv_abund = 1)
  )

  expect_equal(sum(sad1lnorm), n_sim)
  expect_equal(length(sad1lnorm), s_pool)

  sad1bs <- sim_sad(sad_type = "bs", sad_coef = list(S = s_pool, N = n_sim))

  expect_equal(sum(sad1bs), n_sim)
  expect_equal(length(sad1bs), s_pool)

  expect_error(sim_sad(s_pool = 1, n_sim = 0))
})


test_that("sim_sad() - respects seed argument", {
  expect_equal(
    sim_sad(sad_type = "lnorm", s_pool = 10L, n_sim = 500L, seed = 42L),
    sim_sad(sad_type = "lnorm", s_pool = 10L, n_sim = 500L, seed = 42L)
  )
  expect_equal(
    sim_sad(sad_type = "bs", sad_coef = list(S = 4L, N = 10L), seed = 42L),
    sim_sad(sad_type = "bs", sad_coef = list(S = 4L, N = 10L), seed = 42L)
  )
  expect_equal(
    sim_sad(
      sad_type = "gamma",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(shape = 1, scale = 1),
      seed = 42L
    ),
    sim_sad(
      sad_type = "gamma",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(shape = 1, scale = 1),
      seed = 42L
    )
  )
  expect_equal(
    sim_sad(
      sad_type = "geom",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(prob = 0.5),
      seed = 42L
    ),
    sim_sad(
      sad_type = "geom",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(prob = 0.5),
      seed = 42L
    )
  )
  expect_equal(
    sim_sad(sad_type = "ls", sad_coef = list(N = 500L, alpha = 3), seed = 42L),
    sim_sad(sad_type = "ls", sad_coef = list(N = 500L, alpha = 3), seed = 42L)
  )
  expect_equal(
    sim_sad(
      sad_type = "nbinom",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(size = 100L, prob = 0.5),
      seed = 42L
    ),
    sim_sad(
      sad_type = "nbinom",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(size = 100L, prob = 0.5),
      seed = 42L
    )
  )
  expect_equal(
    sim_sad(
      sad_type = "pareto",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(shape = 2),
      seed = 42L
    ),
    sim_sad(
      sad_type = "pareto",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(shape = 2),
      seed = 42L
    )
  )
  expect_equal(
    sim_sad(sad_type = "poilog", s_pool = 10L, n_sim = 500L, seed = 42L),
    sim_sad(sad_type = "poilog", s_pool = 10L, n_sim = 500L, seed = 42L)
  )
  expect_equal(
    sim_sad(
      sad_type = "power",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(s = 2),
      seed = 42L
    ),
    sim_sad(
      sad_type = "power",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(s = 2),
      seed = 42L
    )
  )
  expect_equal(
    sim_sad(
      sad_type = "weibull",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(shape = 2),
      seed = 42L
    ),
    sim_sad(
      sad_type = "weibull",
      s_pool = 10L,
      n_sim = 500L,
      sad_coef = list(shape = 2),
      seed = 42L
    )
  )
})
