## structure ----
test_that("jitter_species() returns a community object", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  result <- jitter_species(simdat, seed = 1L)
  expect_s3_class(result, "community")
  expect_true(all(is.element(
    c("census", "x_min_max", "y_min_max"),
    names(result)
  )))
})

test_that("jitter_species() preserves row count", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  expect_equal(
    nrow(jitter_species(simdat, seed = 1L)$census),
    nrow(simdat$census)
  )
})

## invariants ----
test_that("jitter_species() preserves species column order", {
  # Use an interleaved census so the bug would have been visible
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  # Shuffle rows to interleave species
  set.seed(1L)
  idx <- sample.int(nrow(simdat$census))
  simdat$census <- simdat$census[idx, ]
  rownames(simdat$census) <- NULL

  result <- jitter_species(simdat, sd = 0.1, seed = 1L)
  expect_equal(
    as.character(result$census$species),
    as.character(simdat$census$species)
  )
})

test_that("jitter_species() applies the same x/y shift to all individuals of a species", {
  simdat <- sim_thomas_community(3L, 30L, fix_s_sim = TRUE, seed = 42L)
  result <- jitter_species(simdat, sd = 0.05, seed = 1L)

  for (sp in levels(simdat$census$species)) {
    idx <- simdat$census$species == sp
    if (sum(idx) < 2L) {
      next
    }
    x_diffs <- result$census$x[idx] - simdat$census$x[idx]
    y_diffs <- result$census$y[idx] - simdat$census$y[idx]
    expect_true(
      all(abs(x_diffs - x_diffs[[1L]]) < 1e-12),
      label = paste0("x shift uniform within species ", sp)
    )
    expect_true(
      all(abs(y_diffs - y_diffs[[1L]]) < 1e-12),
      label = paste0("y shift uniform within species ", sp)
    )
  }
})

test_that("jitter_species() coordinates actually change", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  result <- jitter_species(simdat, sd = 0.1, seed = 1L)
  expect_false(isTRUE(all.equal(result$census$x, simdat$census$x)))
})

test_that("jitter_species() with sd = 0 leaves coordinates unchanged", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  result <- jitter_species(simdat, sd = 0, seed = 1L)
  expect_equal(result$census$x, simdat$census$x)
  expect_equal(result$census$y, simdat$census$y)
})

## reproducibility ----
test_that("jitter_species() seed produces identical results", {
  simdat <- sim_thomas_community(5L, 20L, seed = 42L)
  r1 <- jitter_species(simdat, seed = 7L)
  r2 <- jitter_species(simdat, seed = 7L)
  expect_equal(r1$census$x, r2$census$x)
  expect_equal(r1$census$y, r2$census$y)
})

## snapshots ----
test_that("snapshot test with Thoms process", {
  simdat <- sim_thomas_community(10L, 10L, seed = 42L)
  expect_snapshot(jitter_species(simdat, seed = 42L))
})

test_that("snapshot test with Poisson process", {
  simdat <- sim_poisson_community(10L, 10L, seed = 42L)
  expect_snapshot(jitter_species(simdat, seed = 42L))
})

test_that("snapshot test wide jitter", {
  simdat <- sim_thomas_community(10L, 10L, seed = 42L)
  expect_snapshot(jitter_species(simdat, sd = 3, seed = 42L))
})

test_that("snapshot test with singletons only", {
  simdat <- sim_thomas_community(10L, 10L, fix_s_sim = TRUE, seed = 42L)
  expect_snapshot(jitter_species(simdat, seed = 42L))
})
