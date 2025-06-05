# Reproduction functions
test_that("replace_individuals() default snapshot with N specified", {
   simdat <- sim_thomas_coords(11L:100L, seed = 42L)
   expect_snapshot(replace_individuals(simdat, N = 10L, seed = 42L))
})

test_that("replace_individuals() default snapshot with R specified", {
   simdat <- sim_thomas_coords(11L:100L, seed = 42L)
   expect_snapshot(replace_individuals(simdat, R = 0.1, seed = 42L))
})
