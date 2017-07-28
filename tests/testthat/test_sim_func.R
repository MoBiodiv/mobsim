context("Simulation functions")

test_that("classes are correct", {
   sad1 <- sim_sad(s_pool = 100, n_sim = 1000)
   expect_is(sad1, "sad")
   expect_is(sim_poisson_coords(sad1), "community")
   expect_is(sim_poisson_community(100, 1000), "community")
   expect_is(sim_thomas_coords(sad1), "community")
   expect_is(sim_thomas_community(100, 1000), "community")
   expect_is(community(runif(100), runif(100), rep("specA",100)), "community")
   expect_is(community_to_sad(sim_poisson_community(100,1000)), "sad")
})

test_that("species richness and abundance are correct", {
   sad1 <- sim_sad(100,1000, fix_s_sim = T)
   expect_equal(length(sad1), 100)
   expect_equal(sum(sad1), 1000)

   sim1 <- sim_thomas_community(100, 1000, fix_s_sim = T)
   expect_equal(nrow(sim1$census), 1000)
   expect_equal(length(table(sim1$census$species)), 100)
})

