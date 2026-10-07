## invariants ----
test_that("drift_x_species() shifts all x by exactly drift", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  d <- 0.15
  result <- drift_x_species(simdat, drift = d)
  expect_equal(result$census$x, simdat$census$x + d)
})

test_that("drift_y_species() shifts all y by exactly drift", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  d <- -0.3
  result <- drift_y_species(simdat, drift = d)
  expect_equal(result$census$y, simdat$census$y + d)
})

test_that("drift_x_species() does not change y coordinates or species", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  result <- drift_x_species(simdat, drift = 0.1)
  expect_equal(result$census$y, simdat$census$y)
  expect_equal(
    as.character(result$census$species),
    as.character(simdat$census$species)
  )
})

test_that("drift_y_species() does not change x coordinates or species", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  result <- drift_y_species(simdat, drift = 0.1)
  expect_equal(result$census$x, simdat$census$x)
  expect_equal(
    as.character(result$census$species),
    as.character(simdat$census$species)
  )
})

test_that("drift_x_species() preserves community class", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  result <- drift_x_species(simdat, drift = 0.1)
  expect_s3_class(result, "community")
})

test_that("drift_y_species() preserves community class", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  result <- drift_y_species(simdat, drift = 0.1)
  expect_s3_class(result, "community")
})

## edge cases ----
test_that("drift_x_species() with zero drift leaves coordinates unchanged", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  result <- drift_x_species(simdat, drift = 0)
  expect_equal(result$census$x, simdat$census$x)
})

test_that("drift_y_species() with zero drift leaves coordinates unchanged", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  result <- drift_y_species(simdat, drift = 0)
  expect_equal(result$census$y, simdat$census$y)
})

test_that("drift_x_species() errors when drift is not provided", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  expect_error(drift_x_species(simdat))
})

test_that("drift_y_species() errors when drift is not provided", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  expect_error(drift_y_species(simdat))
})

## snapshots ----
test_that("snapshot test with drift", {
  simdat <- sim_thomas_community(10L, 10L, fix_s_sim = TRUE, seed = 42L)
  expect_snapshot(drift_x_species(simdat, drift = 0.1))
  expect_snapshot(drift_y_species(simdat, drift = 0.1))
})
