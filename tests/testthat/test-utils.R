# tests utility functions
## torusify() — invariants ----
test_that("torusify.mother_map() output stays within x_min_max and y_min_max", {
  mpcoords <- create_motherpoints(s_local = 20L, seed = 42L)
  mpcoordsJ <- jitter_motherpoints(mpcoords, sd = 5, seed = 1L)
  result <- torusify(mpcoordsJ)
  all_x <- unlist(result$xmother)
  all_y <- unlist(result$ymother)
  expect_true(all(all_x >= result$x_min_max[[1L]] & all_x < result$x_min_max[[2L]]))
  expect_true(all(all_y >= result$y_min_max[[1L]] & all_y < result$y_min_max[[2L]]))
})

test_that("torusify.community() output stays within x_min_max and y_min_max", {
  simdat <- sim_thomas_community(5L, 50L, seed = 42L)
  simdatJ <- jitter_species(simdat, sd = 5, seed = 1L)
  result <- torusify(simdatJ)
  expect_true(all(
    result$census$x >= result$x_min_max[[1L]] &
    result$census$x < result$x_min_max[[2L]]
  ))
  expect_true(all(
    result$census$y >= result$y_min_max[[1L]] &
    result$census$y < result$y_min_max[[2L]]
  ))
})

test_that("torusify.numeric() output stays within range", {
  vals <- stats::runif(50L, min = -3, max = 3)
  result <- torusify(vals, range = c(0, 1))
  expect_true(all(result >= 0 & result < 1))
})

## torusify.data.frame() — invariants ----
test_that("torusify.data.frame() output stays within x_min_max and y_min_max", {
  set.seed(1L)
  x_min_max <- c(0, 1)
  y_min_max <- c(0, 1)
  df <- data.frame(
    xmother = stats::runif(30L, min = -1, max = 2),
    ymother = stats::runif(30L, min = -1, max = 2)
  )
  result <- torusify(df, x_min_max = x_min_max, y_min_max = y_min_max)
  expect_true(all(result$xmother >= x_min_max[[1L]] & result$xmother < x_min_max[[2L]]))
  expect_true(all(result$ymother >= y_min_max[[1L]] & result$ymother < y_min_max[[2L]]))
})

## torusify.data.frame() — edge cases ----
test_that("torusify.data.frame() errors when xmother column is missing", {
  df <- data.frame(ymother = stats::runif(10L))
  expect_error(
    torusify(df, x_min_max = c(0, 1), y_min_max = c(0, 1)),
    "xmother"
  )
})

test_that("torusify.data.frame() errors when ymother column is missing", {
  df <- data.frame(xmother = stats::runif(10L))
  expect_error(
    torusify(df, x_min_max = c(0, 1), y_min_max = c(0, 1)),
    "ymother"
  )
})

## community2ppp() — structure ----
test_that("community2ppp() output is class ppp", {
  simdat <- sim_thomas_coords(11L:30L, seed = 42L)
  result <- community2ppp(simdat)
  expect_s3_class(result, "ppp")
})

test_that("community2ppp() point count matches community individual count", {
  simdat <- sim_thomas_coords(11L:30L, seed = 42L)
  result <- community2ppp(simdat)
  expect_equal(result$n, nrow(simdat$census))
})

## create_random_ID() — seed and structure ----
test_that("create_random_ID() works in a fresh-seed context (no prior RNG)", {
  # Simulate a fresh session by removing .Random.seed temporarily
  old_seed <- if (exists(".Random.seed")) .Random.seed else NULL
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  rm(".Random.seed", envir = .GlobalEnv)
  expect_no_error(create_random_ID(seed = 1L))
})

test_that("create_random_ID() n > 1 returns correct length and format", {
  result <- create_random_ID(n = 5L, seed = 1L)
  expect_length(result, 5L)
  expect_true(all(nchar(result) == 10L))
  expect_true(all(grepl(pattern = "^[A-Z]{5}[0-9]{4}[A-Z]$", x = result)))
})

## torusify() — snapshots ----
test_that("torusify() - mother_map method - default snapshot 0-1 range", {
  mpcoords <- create_motherpoints(s_local = 20L, seed = 42L)
  mpcoordsJ <- jitter_motherpoints(mpcoords, sd = 3, seed = 42L)
  expect_snapshot(torusify(mpcoordsJ))
})

test_that("torusify() - mother_map method - wild range snapshot", {
  mpcoords <- create_motherpoints(
    s_local = 20L,
    xrange = c(8, 10),
    yrange = c(-1, 1),
    seed = 42L
  )
  mpcoordsJ <- jitter_motherpoints(mpcoords, sd = 3, seed = 42L)
  expect_snapshot(torusify(mpcoordsJ))
})

test_that("torusify() - community method - wild range snapshot", {
  simdat <- sim_thomas_community(
    10L,
    10L,
    xrange = c(8, 10),
    yrange = c(-1, 1),
    seed = 42L
  )
  simdatJ <- jitter_species(simdat, sd = 3, seed = 42L)
  expect_snapshot(torusify(simdatJ))
})

test_that("torusify() - data.frame method - wild range snapshot", {
  set.seed(42L)
  x_min_max <- c(8, 10)
  y_min_max <- c(-1, 1)
  mpcoords <- data.frame(
    xmother = stats::runif(
      20L,
      min = x_min_max[1L] - 1,
      max = x_min_max[2L] + 1
    ),
    ymother = stats::runif(
      20L,
      min = y_min_max[1L] - 1,
      max = y_min_max[2L] + 1
    )
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
  simdat <- sim_thomas_coords(11L:100L)
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

## local_seed() ----
test_that("seeded functions leave the caller's RNG stream untouched", {
  set.seed(seed = 1L)
  expected <- stats::runif(n = 3L)
  set.seed(seed = 1L)
  sim_poisson_community(s_pool = 5L, n_sim = 50L, seed = 2L)
  expect_identical(stats::runif(n = 3L), expected)
})

test_that("seeded functions do not create .Random.seed", {
  old_seed <- get0(x = ".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(
    if (!is.null(old_seed)) {
      assign(x = ".Random.seed", value = old_seed, envir = .GlobalEnv)
    },
    add = TRUE
  )
  rm(list = ".Random.seed", envir = .GlobalEnv)
  jitter_species(sim_poisson_coords(abund_vec = 1L:5L, seed = 1L), sd = 0.1, seed = 2L)
  expect_false(exists(x = ".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
