# Reproduction functions

## edge cases ----
test_that("replace_individuals() errors when both N and R are NULL", {
  simdat <- sim_thomas_coords(11L:30L, seed = 42L)
  expect_error(
    replace_individuals(simdat),
    regexp = "N.*R|R.*N",
    ignore.case = TRUE
  )
})

## structure ----
test_that("replace_individuals() output is class community", {
  simdat <- sim_thomas_coords(11L:30L, seed = 42L)
  result <- replace_individuals(simdat, N = 5L, seed = 42L)
  expect_s3_class(result, "community")
})

test_that("replace_individuals() preserves individual count", {
  simdat <- sim_thomas_coords(11L:30L, seed = 42L)
  result <- replace_individuals(simdat, N = 5L, seed = 42L)
  expect_equal(nrow(result$census), nrow(simdat$census))
})

test_that("replace_individuals() species column is a factor", {
  simdat <- sim_thomas_coords(11L:30L, seed = 42L)
  result <- replace_individuals(simdat, N = 5L, seed = 42L)
  expect_true(is.factor(result$census$species))
})

test_that("replace_individuals() x and y coordinates are unchanged", {
  simdat <- sim_thomas_coords(11L:30L, seed = 42L)
  result <- replace_individuals(simdat, N = 5L, seed = 42L)
  expect_equal(result$census$x, simdat$census$x)
  expect_equal(result$census$y, simdat$census$y)
})

## invariants ----
test_that("replace_individuals() only assigns existing species", {
  simdat <- sim_thomas_coords(11L:30L, seed = 42L)
  result <- replace_individuals(simdat, N = 5L, seed = 42L)
  original_species <- levels(simdat$census$species)
  result_species <- levels(result$census$species)
  expect_true(all(is.element(result_species, original_species)))
})

## reproducibility ----
test_that("replace_individuals() seed produces identical results", {
  simdat <- sim_thomas_coords(11L:30L, seed = 42L)
  r1 <- replace_individuals(simdat, N = 5L, seed = 7L)
  r2 <- replace_individuals(simdat, N = 5L, seed = 7L)
  expect_equal(as.character(r1$census$species), as.character(r2$census$species))
})

## snapshots ----
test_that("replace_individuals() default snapshot with N specified", {
  simdat <- sim_thomas_coords(11L:100L, seed = 42L)
  expect_snapshot(replace_individuals(simdat, N = 10L, seed = 42L))
})

test_that("replace_individuals() default snapshot with R specified", {
  simdat <- sim_thomas_coords(11L:100L, seed = 42L)
  expect_snapshot(replace_individuals(simdat, R = 0.1, seed = 42L))
})
