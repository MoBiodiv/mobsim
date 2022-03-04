test_that("sim_poisson_community() - default argument values are as expected", {
   expect_snapshot_output(sim_poisson_community(s_pool = 3L, n_sim = 100L, seed = 42L))
})

test_that("sim_poisson_coords() - default argument values are as expected", {
   expect_snapshot_output(sim_poisson_coords(abund_vec = 1:4, seed = 42L))
})

test_that("sim_poisson_community() calls sim_sad() correctly", {
   s_pool <- 10L
   n_sim <- 100L

   sad1 <- sim_sad(s_pool, n_sim, fix_s_sim = TRUE, seed = 42L)
   sim1 <- sim_poisson_community(s_pool, n_sim, fix_s_sim = TRUE, seed = 42L)

   expect_equal(sum(sad1), nrow(sim1$census))
   expect_equal(length(sad1), length(table(sim1$census$species)))
})
