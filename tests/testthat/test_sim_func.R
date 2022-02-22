# Testing simulation functions

test_that("classes are correct", {
   sad1 <- sim_sad(s_pool = 3L, n_sim = 10L)
   expect_is(sad1, "sad")
   expect_is(sad1, "integer")

   xmother <- lapply(1:length(sad1), function(x) runif(1, 0, 1))
   ymother <- lapply(1:length(sad1), function(x) runif(1, 0, 1))

   expect_is(sim_poisson_coords(sad1), "community")
   expect_is(sim_poisson_community(3L, 10L), "community")
   expect_is(sim_thomas_coords(sad1), "community")
   expect_is(sim_thomas_coords(sad1, mother_points = 0), "community")
   expect_is(sim_thomas_coords(sad1, mother_points = 100), "community")
   expect_is(sim_thomas_coords(sad1, mother_points = 4.2), "community")
   expect_is(sim_thomas_coords(sad1, xmother = xmother, ymother = ymother), "community")
   expect_is(sim_thomas_coords(sad1), "community")
   expect_is(sim_thomas_community(3L, 10L), "community")
   expect_is(community(runif(10L), runif(10L), rep("specA", 10L)), "community")
   expect_is(community_to_sad(sim_poisson_community(3L, 10L)), "sad")
})

test_that("sim_sad() - correct assertions", {
   expect_error(sim_sad(n_sim = 100L, sad_type = "lnorm"))
   expect_error(sim_sad(s_pool = NA, n_sim = 100L, sad_type = "lnorm"), "The argument s_pool is mandatory for the selected sad and has to be a positive integer number.")
   expect_error(sim_sad(s_pool = 3.5, n_sim = 100L, sad_type = "lnorm"), "s_pool has to be a positive integer number")
   expect_error(sim_sad(s_pool = 0L, n_sim = 100L, sad_type = "lnorm"), "s_pool has to be a positive integer number")

   expect_error(sim_sad(s_pool = 3L, sad_type = "lnorm"))
   expect_error(sim_sad(s_pool = 3L, n_sim = NA, sad_type = "lnorm"))
   expect_error(sim_sad(s_pool = 3L, n_sim = 100.5, sad_type = "lnorm"), "n_sim has to be a positive integer number")
   expect_error(sim_sad(s_pool = 3L, n_sim = 0L, sad_type = "lnorm"), "n_sim has to be a positive integer number")

   expect_error(sim_sad(s_pool = 3L, n_sim = 100L, sad_type = "lnorm", sad_coef = 0.1), "coef must be a named list!")
   expect_error(sim_sad(s_pool = 3L, n_sim = 100L, sad_type = "wrong_sad_type"))

   expect_error(sim_sad(sad_type = "bs", sad_coef = list(S = 4L)))

   expect_warning(sim_sad(s_pool = 3L, sad_type = "bs", sad_coef = list(S = 4L, N = 10L)))
   expect_warning(sim_sad(n_sim = 10L, sad_type = "bs", sad_coef = list(S = 4L, N = 10L)))
})


test_that("sim_sad() - results are as expected", {
   expect_equal(names(sim_sad(s_pool = 3L, n_sim = 10L))[1L], "species_1")
   n_sim <- 10L
   s_pool <- 3L

   sadNorm <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = "lnorm", fix_s_sim = TRUE)
   expect_equal(sum(sadNorm), n_sim)
   expect_equal(length(sadNorm), s_pool)

   sadBs <- sim_sad(sad_type = "bs", sad_coef = list(S = s_pool, N = n_sim), fix_s_sim = TRUE)
   expect_equal(sum(sadBs), n_sim)
   expect_equal(length(sadBs), s_pool)

   sadGamma <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = "gamma", sad_coef = list(shape = 1, scale = 1), fix_s_sim = TRUE)
   expect_equal(sum(sadGamma), n_sim)
   expect_equal(length(sadGamma), s_pool)

   sadGeom <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = "geom",  sad_coef = list(prob = 0.5), fix_s_sim = TRUE)
   expect_equal(sum(sadGeom), n_sim)
   expect_equal(length(sadGeom), s_pool)

   sadLs <- sim_sad(sad_type = "ls", sad_coef = list(N = n_sim, alpha = 3), fix_s_sim = TRUE)
   expect_equal(sum(sadLs), n_sim)

   # expect_equal(sum(sim_sad(sad_type = "mzsm", sad_coef = list(n = 3L, J = n_sim, theta = 4))), n_sim)

   sadNbinom <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = "nbinom", sad_coef = list(size = 100L, prob = 0.5), fix_s_sim = TRUE)
   expect_equal(sum(sadNbinom), n_sim)
   expect_equal(length(sadNbinom), s_pool)

   sadPareto <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = "pareto", sad_coef = list(shape = 2), fix_s_sim = TRUE)
   expect_equal(sum(sadPareto), n_sim)
   expect_equal(length(sadPareto), s_pool)

   sadPoilog <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = "poilog", fix_s_sim = TRUE)
   expect_equal(sum(sadPoilog), n_sim)
   expect_equal(length(sadPoilog), s_pool)

   sadPower <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = "power", sad_coef = list(s = 2), fix_s_sim = TRUE)
   expect_equal(sum(sadPower), n_sim)
   expect_equal(length(sadPower), s_pool)

   # expect_equal(sum(sim_sad(s_pool = 3L, n_sim = n_sim, sad_type = "powbend", sad_coef = list(s = 2, omega = 3))), n_sim)

   sadWeibull <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = "weibull", sad_coef = list(shape = 2), fix_s_sim = TRUE)
   expect_equal(sum(sadWeibull), n_sim)
   expect_equal(length(sadWeibull), s_pool)
})

test_that("sim_sad() - edge case s_pool = 1", {
   n_sim <- 100L
   s_pool <- 1L
   sad1lnorm <- sim_sad(s_pool = s_pool, n_sim = n_sim, sad_type = "lnorm", sad_coef = list(cv_abund = 1))

   expect_equal(sum(sad1lnorm), n_sim)
   expect_equal(length(sad1lnorm), s_pool)

   sad1bs <- sim_sad(sad_type = "bs", sad_coef = list(S = s_pool, N = n_sim))

   expect_equal(sum(sad1bs), n_sim)
   expect_equal(length(sad1bs), s_pool)

   expect_error(sim_sad(s_pool = 1, n_sim = 0))
})


test_that("sim_sad() - respects seed argument", {
   expect_equal(
      sim_sad(sad_type = "lnorm", s_pool = 10L, n_sim = 500L, seed = 42L),
      sim_sad(sad_type = "lnorm", s_pool = 10L, n_sim = 500L, seed = 42L)
   )
   expect_equal(
      sim_sad(sad_type = "bs", sad_coef = list(S = 4L, N = 10L), seed = 42L),
      sim_sad(sad_type = "bs", sad_coef = list(S = 4L, N = 10L), seed = 42L)
   )
   expect_equal(
      sim_sad(sad_type = "gamma", s_pool = 10L, n_sim = 500L, sad_coef = list(shape = 1, scale = 1), seed = 42L),
      sim_sad(sad_type = "gamma", s_pool = 10L, n_sim = 500L, sad_coef = list(shape = 1, scale = 1), seed = 42L)
   )
   expect_equal(
      sim_sad(sad_type = "geom",  s_pool = 10L, n_sim = 500L, sad_coef = list(prob = 0.5), seed = 42L),
      sim_sad(sad_type = "geom",  s_pool = 10L, n_sim = 500L, sad_coef = list(prob = 0.5), seed = 42L)
   )
   expect_equal(
      sim_sad(sad_type = "ls", sad_coef = list(N = 500L, alpha = 3), seed = 42L),
      sim_sad(sad_type = "ls", sad_coef = list(N = 500L, alpha = 3), seed = 42L)
   )
   expect_equal(
      sim_sad(sad_type = "nbinom", s_pool = 10L, n_sim = 500L, sad_coef = list(size = 100L, prob = 0.5), seed = 42L),
      sim_sad(sad_type = "nbinom", s_pool = 10L, n_sim = 500L, sad_coef = list(size = 100L, prob = 0.5), seed = 42L)
   )
   expect_equal(
      sim_sad(sad_type = "pareto", s_pool = 10L, n_sim = 500L, sad_coef = list(shape = 2), seed = 42L),
      sim_sad(sad_type = "pareto", s_pool = 10L, n_sim = 500L, sad_coef = list(shape = 2), seed = 42L)
   )
   expect_equal(
      sim_sad(sad_type = "poilog", s_pool = 10L, n_sim = 500L, seed = 42L),
      sim_sad(sad_type = "poilog", s_pool = 10L, n_sim = 500L, seed = 42L)
   )
   expect_equal(
      sim_sad(sad_type = "power", s_pool = 10L, n_sim = 500L, sad_coef = list(s = 2), seed = 42L),
      sim_sad(sad_type = "power", s_pool = 10L, n_sim = 500L, sad_coef = list(s = 2), seed = 42L)
   )
   expect_equal(
      sim_sad(sad_type = "weibull", s_pool = 10L, n_sim = 500L, sad_coef = list(shape = 2), seed = 42L),
      sim_sad(sad_type = "weibull", s_pool = 10L, n_sim = 500L, sad_coef = list(shape = 2), seed = 42L)
   )
})

test_that("sim_thomas_community() calls sim_sad() correctly", {
   s_pool <- 10L
   n_sim <- 100L

   sad1 <- sim_sad(s_pool, n_sim, fix_s_sim = TRUE, seed = 42L)
   sim1 <- sim_thomas_community(s_pool, n_sim, fix_s_sim = TRUE, seed = 42L)

   expect_equal(sum(sad1), nrow(sim1$census))
   expect_equal(length(sad1), length(table(sim1$census$species)))
})

test_that("sim_thomas_coords() - handles wrong mother_points parametres", {
   expect_error(sim_thomas_coords(sad1, mother_points = -1))  # mother_points is negative
   expect_error(sim_thomas_coords(sad1, mother_points = sample(size = 100, -2:2, replace = TRUE)))  # mother_points countains negative values
   expect_error(sim_thomas_coords(sad1, mother_points = rep(1, 4)))  # mother_points is too short
   expect_error(sim_thomas_coords(sad1, mother_points = rep(1, 1000)))  # mother_points is too long
})

test_that("sim_thomas_coords() - handles wrong xmother and ymother parametres", {
   sad2 <- sim_sad(s_pool = 2, n_sim = 100)
   xmother <- lapply(1:length(sad2), function(x) runif(2, 0, 1))

   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother = list(c(0.2, NA), c(0.3, NA))))
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother = list(NA, c(0.2, 0.3))))
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother = list(2, c(0.2, 0.3))))
})

test_that("sim_thomas_coords() - handles wrong xrange and yrange parametres", {
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother = xmother, xrange = c(0, 0.1)))   # xmother and ymother outside of range
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother = xmother, xrange = data.frame(c(0,1), c(0,0.1))))   # xrange and y range have different class
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother = xmother,
                                  xrange = data.frame(c(0,1,2), c(0,1,2)),
                                  yrange = data.frame(c(0,1,2), c(0,1,2))))  # xrange and yrange have too many columns
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother = xmother,
                                  xrange = data.frame(c(0,1), c(0,1), c(0,1)),
                                  yrange = data.frame(c(0,1), c(0,1), c(0,1))))  # xrange and yrange have too many rows
})

test_that("sim_thomas_coords() - throws warnings when expected", {
   sad1 <- sim_sad(s_pool = 3L, n_sim = 10L)
   xmother <- lapply(1L:length(sad1), function(x) runif(1, 0, 1))
   ymother <- lapply(1L:length(sad1), function(x) runif(1, 0, 1))
   expect_warning(sim_thomas_coords(sad1, mother_points = 1, xmother = xmother, ymother = ymother))
   expect_warning(sim_thomas_coords(sad1, mother_points = 1, cluster_points = 2))
})

test_that("rThomas_rcpp() - behaves as expected", {
   expect_is(rThomas_rcpp(20, n_mother_points = 1, xmother = .5, ymother = .5, sigma = 0.2), "data.frame")
   expect_equal(nrow(rThomas_rcpp(20, n_mother_points = 1, xmother = .5, ymother = .5, sigma = 0.2)), 20)
   expect_error(rThomas_rcpp(20, n_mother_points = 2, sigma = 0.2))
})


if (FALSE) {
   sad2 <- sim_sad(4, 100)
   mother_points <-  c(1,-1,1,1)
   plot(sim_thomas_coords(sad2, mother_points = mother_points))

   xmother <- lapply(1:length(sad2), function(x) runif(2, 0, 1))
   ymother <- lapply(1:length(sad2), function(x) runif(2, 0, 1))
   sim_thomas_coords(sad2, xmother = xmother, ymother = ymother, xrange = c(0, 0.1), yrange = c(0, 0.1))
   # cluster_points
   # NA, class, length

   # xrange yrange
   ## solved
   # if xrange and yrange are not of the same class. NOW error
   # if xrange and yrange are not of the good dimension. NOW error
   # if xrange and yrange are

   # xmother ymother
   ## solved
   # if xmother or ymother coordinates are outside of xrange or yrange, R crashes. NOW error
   # if coordinates are given but mixed with NAs, coordinates are recycled.

   # mother_points
   # if mother_points is numeric, it is truncated to integer by rThomas_rcpp
   # if mother_points > values of abund_vec, each mother point has 1 or 0 individual and the distribution is random
   ## Solved
   # if mother_points is of the same length as abund_vec and one or more value < 0, no clustering for this one but correct clustering for the others. No warning() -> NOW error
   # if mother_points is length 1 and value < 0, no clustering. No warning() -> NOW error
   # if mother_points countains one NA, other values are not taken into account and the function computes mother_points per species on its own. No warning. -> NOW error
   # if length(mother_points) < length(abund_vec), values are not taken into account and the function computes mother_points per species on its own. No warning -> NOW error
   # if length(mother_points) > length(abund_vec), only values 1 to length(abund_vec) are used. -> NOW error
}
