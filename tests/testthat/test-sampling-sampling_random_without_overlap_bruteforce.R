test_that("warns if there is overlap", {
   expect_warning(
      sampling_random_bruteforce(n_quadrats = 10L,
                                 xmin = 0, xmax = 1,
                                 ymin = 0, ymax = 1,
                                 min_dist = 1
      )
   )
})

test_that("classes are correct", {
   expect_type(
      sampling_random_bruteforce(n_quadrats = 2L,
                                 xmin = 0, xmax = 1,
                                 ymin = 0, ymax = 1,
                                 min_dist = sqrt(2 * 0.01)
      ), "list")
   expect_s3_class(
      sampling_random_bruteforce(n_quadrats = 2L,
                                 xmin = 0, xmax = 1,
                                 ymin = 0, ymax = 1,
                                 min_dist = sqrt(2 * 0.01)
      ), "data.frame")
})

test_that("dimensions are correct", {
   n_quadrats = 5L
   xy_dat <- sampling_random_bruteforce(n_quadrats = n_quadrats,
                                        xmin = 0, xmax = 1,
                                        ymin = 0, ymax = 1,
                                        min_dist = sqrt(2 * 0.01)
   )

   expect_equal(dim(xy_dat), c(n_quadrats, 2L))
})

test_that("samples are within range", {
   skip("sample_quadrats() ensures that this does not happen")
   xmin = 0
   xmax = 1
   ymin = 0
   ymax = 1
   quadrat_size =  2

   xy_dat <- sampling_random_bruteforce(n_quadrats = 2L,
                                        xmin = xmin, xmax = xmax,
                                        ymin = ymin, ymax = ymax,
                                        min_dist = sqrt(2 * quadrat_size^2)
   )

   expect_gte(min(xy_dat$x), xmin)
   expect_lte(max(xy_dat$x) + quadrat_size, xmax)
   expect_gte(min(xy_dat$y), ymin)
   expect_lte(max(xy_dat$y) + quadrat_size, ymax)
})


test_that("seed parameter is respected", {
   sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)

   expect_equal(
      sampling_random_bruteforce(n_quadrats = 5L,
                                 xmin = 0, xmax = 1,
                                 ymin = 0, ymax = 1,
                                 min_dist = sqrt(2 * 0.01),
                                 seed = 42L),
      sampling_random_bruteforce(n_quadrats = 5L,
                                 xmin = 0, xmax = 1,
                                 ymin = 0, ymax = 1,
                                 min_dist = sqrt(2 * 0.01),
                                 seed = 42L)
   )
})
