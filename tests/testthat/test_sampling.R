# Tests for sample_quadrats are encapsulated inside the function testing_sampling()
# because tests are repeated to check the 3 methods of placing quadrats.

testing_sampling <- function(method = "random", avoid_overlap = TRUE, mock_spatstat_absence = FALSE) {

   if (method %in% c("grid","transect")) {
      context(paste("Testing sample_quadrats() -", method))
   } else {
      if (!avoid_overlap) {
         context("Testing sample_quadrats() - random with overlap")
      } else {
         if (mock_spatstat_absence) {
            mockery::stub(where = sample_quadrats, what = "requireNamespace", how = FALSE, depth = 1L)
            context("Testing sample_quadrats() - random without overlap - runif")
         } else {
            context("Testing sample_quadrats() - random without overlap - spatstat")
         }
      }
   }
   test_that("fails correctly", {
      sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)

      expect_error(sample_quadrats(data.frame(c(1:5), c(1:5)), method = method, avoid_overlap = avoid_overlap))
      expect_error(sample_quadrats(sim_com1, quadrat_area = 1, method = method, avoid_overlap = avoid_overlap))
      expect_error(sample_quadrats(sim_com1, quadrat_area = 0.01, method = "foo", avoid_overlap = avoid_overlap))
      if (method == "transect") {
         expect_error(sample_quadrats(sim_com1, quadrat_area = 0.9, n_quadrats = 2L,
                                      plot = FALSE, method = method, avoid_overlap = avoid_overlap))
      } else {
         expect_warning(sample_quadrats(sim_com1, quadrat_area = 0.9, n_quadrats = 2L,
                                        plot = FALSE, method = method, avoid_overlap = avoid_overlap))
      }
   })
   test_that("sampling_transects() - fails correctly", {
      expect_error(
         sampling_transects(n_quadrats = 2L,
                            xmin = 0, xmax = 1,
                            ymin = 0, ymax = 1,
                            x0 = 0, y0 = 0, delta_x = 1, delta_y = 0.1,
                            quadrat_size = sqrt(0.01)
         )
      )

      expect_error(
         sampling_transects(n_quadrats = 2L,
                            xmin = 0, xmax = 1,
                            ymin = 0, ymax = 1,
                            x0 = 0, y0 = 0, delta_x = 0.1, delta_y = 1,
                            quadrat_size = sqrt(0.01)
         )
      )

      expect_warning(
         sampling_transects(n_quadrats = 10L,
                            xmin = 0, xmax = 1,
                            ymin = 0, ymax = 1,
                            x0 = 0.4, y0 = 0.4, delta_x = 0.01, delta_y = 0.01,
                            quadrat_size = 0.1
         )
      )
   })

   test_that("sampling_grids() - fails correctly", {
      expect_error(
         sampling_grids(n_quadrats = 2L,
                        xmin = 0.1, xmax = 1, # xmin is larger than x0
                        ymin = 0, ymax = 1,
                        x0 = 0, y0 = 0, delta_x = 0.1, delta_y = 0.1,
                        quadrat_size = 0.01
         )
      )

      expect_error(
         sampling_grids(n_quadrats = 2L,
                        xmin = 0, xmax = 1,
                        ymin = 0.1, ymax = 1, # ymin is larger than y0
                        x0 = 0, y0 = 0, delta_x = 0.1, delta_y = 0.1,
                        quadrat_size = 0.01
         )
      )

      expect_warning(
         sampling_grids(n_quadrats = 10L,
                        xmin = 0, xmax = 1,
                        ymin = 0, ymax = 1,
                        x0 = 0.4, y0 = 0.4, delta_x = 0.01, delta_y = 0.01,
                        quadrat_size = 0.1
         )
      )
   })


   test_that("classes are correct", {
      sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
      comm_mat1 <- sample_quadrats(sim_com1, plot = FALSE, n_quadrats = 2L,
                                   method = method, avoid_overlap = avoid_overlap)

      expect_type(comm_mat1, "list")
      expect_s3_class(comm_mat1[[1]], "data.frame")
      expect_s3_class(comm_mat1[[2]], "data.frame")
   })

   test_that("dimensions are correct", {
      S = 5L
      sim_com1 <- sim_poisson_community(s_pool = S, n_sim = 50L, fix_s_sim = TRUE)
      n_quadrats <- 10L
      comm_mat <- sample_quadrats(
         sim_com1, n_quadrats = n_quadrats,
         method = method, plot = FALSE,
         avoid_overlap = avoid_overlap
      )

      expect_vector(comm_mat, size = 2L)
      expect_equal(dim(comm_mat[[1]]), c(n_quadrats, S))
      expect_equal(dim(comm_mat[[2]]), c(n_quadrats, 2L))

   })

   test_that("samples are within range", {
      quadrat_area <- 0.002
      sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
      comm_mat1 <- sample_quadrats(
         sim_com1, n_quadrats = 10L, quadrat_area = quadrat_area,
         method = method, plot = FALSE,
         avoid_overlap = avoid_overlap
      )

      expect_true(sim_com1$x_min_max[1] <=  min(comm_mat1[[2]][,"x"]) && max(comm_mat1[[2]][,"x"]) + sqrt(quadrat_area) <= sim_com1$x_min_max[2])
      expect_true(sim_com1$y_min_max[1] <=  min(comm_mat1[[2]][,"y"]) && max(comm_mat1[[2]][,"y"]) + sqrt(quadrat_area) <= sim_com1$y_min_max[2])
   })

   test_that("seed parameter is respected", {
      sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)

      expect_equal(
         sample_quadrats(sim_com1, seed = 42L, plot = FALSE, n_quadrats = 5L, method = method, avoid_overlap = avoid_overlap),
         sample_quadrats(sim_com1, seed = 42L, plot = FALSE, n_quadrats = 5L, method = method, avoid_overlap = avoid_overlap)
      )
   })

   test_that("sample_quadrats() - edge case where n_quadrats == 1", {
      sim_com1 <- sim_poisson_community(s_pool = 5L, n_sim = 50L)
      quadrat_area <- 0.01
      expect_equal(
         as.numeric(
            sample_quadrats(sim_com1, seed = 42L, plot = FALSE, n_quadrats = 1L, quadrat_area = quadrat_area,
                            method = method, avoid_overlap = avoid_overlap)$xy_dat
         ),
         as.numeric(
            sampling_one_quadrat(xmin = sim_com1$x_min_max[1L], xmax = sim_com1$x_min_max[2L] - sqrt(quadrat_area),
                                 ymin = sim_com1$y_min_max[1L], ymax = sim_com1$y_min_max[2L] - sqrt(quadrat_area), seed = 42L)
         )
      )
   })


} # end of testing_sampling()


# if avoid_overlap = FALSE, random sampling is used
testing_sampling(method = "random", avoid_overlap = FALSE)

# if avoid_overlap = TRUE and package spatstat.random is installed (hence requireNamespace = TRUE inside sample_quadrats), spatstat.random is always used BUT we also want to test without the package so a stub is created to force requireNamespace("spatstat.random", quietly = TRUE) to return FALSE.
if (requireNamespace("spatstat.random", quietly = TRUE)) {

   testing_sampling(method = "random", avoid_overlap = TRUE)
   testing_sampling(method = "random", avoid_overlap = TRUE, mock_spatstat_absence = TRUE)

} else {

   testing_sampling(method = "random", avoid_overlap = TRUE)

}

# If method is not "random", avoid_overlap is ignored.
testing_sampling(method = "transect")

testing_sampling(method = "grid")

