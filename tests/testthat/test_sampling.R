context("Testing sample_quadrats general behaviour")

test_that("fails correctly", {
   expect_error(
      sample_quadrats(
         data.frame(c(1:5), c(1:5)), n_quadrats = 10, quadrat_area = 0.001,
         method = "grid", plot = FALSE
      ))
})

test_that("classes are correct", {
   S <-  5L
   N <-  50L
   n_quadrats <- 10L
   quadrat_area <- 0.002
   sim_com1 <- sim_poisson_community(s_pool = S, n_sim = N)
   comm_mat1 <- sample_quadrats(
      sim_com1, n_quadrats = n_quadrats, quadrat_area = quadrat_area,
      method = "grid", plot = FALSE
   )

   expect_is(comm_mat1, "list")
   expect_is(comm_mat1[[1]], "data.frame")
   expect_is(comm_mat1[[2]], "data.frame")
})

test_that("dimensions are correct", {
   S <-  5L
   N <-  50L
   n_quadrats <- 10L
   quadrat_area <- 0.002
   sim_com1 <- sim_poisson_community(s_pool = S, n_sim = N)
   comm_mat1 <- sample_quadrats(
      sim_com1, n_quadrats = n_quadrats, quadrat_area = quadrat_area,
      method = "grid", plot = FALSE
   )

   expect_vector(comm_mat1, size = 2L)
   expect_equal(dim(comm_mat1[[1]]), c(n_quadrats, S))
   expect_equal(dim(comm_mat1[[2]]), c(n_quadrats, 2))
})

test_that("samples are within range", {
   S <-  5L
   N <-  50L
   n_quadrats <- 10L
   quadrat_area <- 0.002
   sim_com1 <- sim_poisson_community(s_pool = S, n_sim = N)
   comm_mat1 <- sample_quadrats(
      sim_com1, n_quadrats = n_quadrats, quadrat_area = quadrat_area,
      method = "grid", plot = FALSE
   )

   expect_true(sim_com1$x_min_max[1] <=  min(comm_mat1[[2]][,"x"]) && max(comm_mat1[[2]][,"x"]) + sqrt(quadrat_area) <= sim_com1$x_min_max[2])
   expect_true(sim_com1$y_min_max[1] <=  min(comm_mat1[[2]][,"y"]) && max(comm_mat1[[2]][,"y"]) + sqrt(quadrat_area) <= sim_com1$y_min_max[2])
})


context("Testing sample_quadrats with spatstat.random installed")
context("Testing sample_quadrats without spatstat.random")

#testing has to happen in 2 cases: with and without spatstat
