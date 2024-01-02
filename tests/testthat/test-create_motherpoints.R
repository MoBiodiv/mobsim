test_that("snapshot test with clustering", {
   expect_snapshot(create_motherpoints(s_local = 100L, seed = 42L))
})

test_that("snapshot test with clustering - no range", {
   expect_snapshot(create_motherpoints(s_local = 100L,
                                       n_motherpoints_range = 1L,
                                       seed = 42L))
})

test_that("snapshot test without clustering", {
   expect_snapshot(create_motherpoints(s_local = 100L,
                                       n_motherpoints_range = 0L:1L,
                                       seed = 42L))
})

test_that("snapshot test plotting method for mother_map class", {
   mpcoords <- create_motherpoints(s_local = 100L,
                                   n_motherpoints_range = 0L:1L,
                                   seed = 42L)
   expect_snapshot_output(plot(mpcoords))
})
