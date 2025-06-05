test_that("snapshot test with drift", {
  simdat <- sim_thomas_community(10L, 10L, fix_s_sim = TRUE, seed = 42L)
  expect_snapshot(drift_x_species(simdat, drift = 0.1))
  expect_snapshot(drift_y_species(simdat, drift = 0.1))
})
