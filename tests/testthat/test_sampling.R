context("Testing sample_quadrats general behaviour")

test_that("fails correctly", {
   sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)

   expect_error(sample_quadrats(data.frame(c(1:5), c(1:5))))
   expect_error(sample_quadrats(sim_com1, quadrat_area = 1))
   expect_error(sample_quadrats(sim_com1, quadrat_area = 0.01, method = "foo"))

   expect_warning(sample_quadrats(sim_com1, quadrat_area = 0.9, n_quadrats = 2L, plot = FALSE))
})

test_that("classes are correct", {
   sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
   comm_mat1 <- sample_quadrats(sim_com1, plot = FALSE, n_quadrats = 2L)

   expect_is(comm_mat1, "list")
   expect_is(comm_mat1[[1]], "data.frame")
   expect_is(comm_mat1[[2]], "data.frame")
})

test_that("dimensions are correct", {
   sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
   n_quadrats <- 10L
   comm_matGrid <- sample_quadrats(
      sim_com1, n_quadrats = n_quadrats,
      method = "grid", plot = FALSE
   )
   comm_matTransect <- sample_quadrats(
      sim_com1, n_quadrats = n_quadrats,
      method = "transect", plot = FALSE
   )
   comm_matRandom <- sample_quadrats(
      sim_com1, n_quadrats = n_quadrats,
      method = "random", plot = FALSE, seed = 2
   )


   expect_vector(comm_matGrid, size = 2L)
   expect_equal(dim(comm_matGrid[[1]]), c(n_quadrats, 5L))
   expect_equal(dim(comm_matGrid[[2]]), c(n_quadrats, 2))

   expect_vector(comm_matTransect, size = 2L)
   expect_equal(dim(comm_matTransect[[1]]), c(n_quadrats, 5L))
   expect_equal(dim(comm_matTransect[[2]]), c(n_quadrats, 2))

   expect_vector(comm_matTransect, size = 2L)
   expect_equal(dim(comm_matTransect[[1]]), c(n_quadrats, 5L))
   expect_equal(dim(comm_matTransect[[2]]), c(n_quadrats, 2))
})

test_that("samples are within range", {
   quadrat_area <- 0.002
   sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
   comm_mat1 <- sample_quadrats(
      sim_com1, n_quadrats = 50L, quadrat_area = quadrat_area,
      method = "grid", plot = FALSE
   )

   expect_true(sim_com1$x_min_max[1] <=  min(comm_mat1[[2]][,"x"]) && max(comm_mat1[[2]][,"x"]) + sqrt(quadrat_area) <= sim_com1$x_min_max[2])
   expect_true(sim_com1$y_min_max[1] <=  min(comm_mat1[[2]][,"y"]) && max(comm_mat1[[2]][,"y"]) + sqrt(quadrat_area) <= sim_com1$y_min_max[2])
})

test_that("seed parameter is respected", {
   sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)

   expect_equal(
      sample_quadrats(sim_com1, seed = 42L, plot = FALSE, n_quadrats = 5L),
      sample_quadrats(sim_com1, seed = 42L, plot = FALSE, n_quadrats = 5L)
   )
})


context("Testing sample_quadrats with spatstat.random installed")
context("Testing sample_quadrats without spatstat.random")
