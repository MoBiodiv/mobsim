context("Simulation functions")

test_that("classes are correct", {
   sad1 <- sim_sad(s_pool = 100, n_sim = 1000)
   xmother <- lapply(1:length(sad1), function(x) runif(1, 0, 1))
   ymother <- lapply(1:length(sad1), function(x) runif(1, 0, 1))
   expect_is(sad1, "sad")
   expect_is(sim_poisson_coords(sad1), "community")
   expect_is(sim_poisson_community(100, 1000), "community")
   expect_is(sim_thomas_coords(sad1), "community")
   expect_is(sim_thomas_coords(sad1, mother_points = 0), "community")
   expect_is(sim_thomas_coords(sad1, mother_points = 100), "community")
   expect_is(sim_thomas_coords(sad1, mother_points = 4.2), "community")
   expect_is(sim_thomas_coords(sad1, mother_points = 1, cluster_points = 2), "community")
   expect_is(sim_thomas_coords(sad1, mother_points = 1, xmother=xmother, ymother=ymother), "community")
   expect_is(sim_thomas_coords(sad1, xmother=xmother, ymother=ymother), "community")
   expect_is(sim_thomas_coords(sad1), "community")
   expect_is(sim_thomas_community(100, 1000), "community")
   expect_is(community(runif(100), runif(100), rep("specA",100)), "community")
   expect_is(community_to_sad(sim_poisson_community(100,1000)), "sad")
})

test_that("The function handles wrong mother_points parametres", {
   expect_error(sim_thomas_coords(sad1, mother_points = -1))  # mother_points is negative
   expect_error(sim_thomas_coords(sad1, mother_points = sample(size=100, -2:2, replace=T)))  # mother_points countains negative values
   expect_error(sim_thomas_coords(sad1, mother_points = rep(1, 4)))  # mother_points is too short
   expect_error(sim_thomas_coords(sad1, mother_points = rep(1, 1000)))  # mother_points is too long
})

test_that("The function handles wrong xmother and ymother parametres", {
   sad2 <- sim_sad(s_pool = 2, n_sim = 100)
   xmother <- lapply(1:length(sad2), function(x) runif(2, 0, 1))
   
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother=list(c(0.2, NA), c(0.3, NA))))
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother=list(NA, c(0.2, 0.3))))
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother=list(2, c(0.2, 0.3))))
})

test_that("The function rThomas_rcpp behaves as expected", {
   expect_is(rThomas_rcpp(20, n_mother_points = 1, xmother=.5, ymother=.5, sigma=0.2), "data.frame")
   expect_equal(nrow(rThomas_rcpp(20, n_mother_points = 1, xmother=.5, ymother=.5, sigma=0.2)), 20)
   expect_error(rThomas_rcpp(20, n_mother_points = 2, sigma=0.2))
})

test_that("The function handles wrong xrange and yrange parametres", {
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother=xmother, xrange=c(0, 0.1)))   # xmother and ymother outside of range
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother=xmother, xrange=data.frame(c(0,1), c(0,0.1))))   # xrange and y range have different class
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother=xmother,
                                  xrange=data.frame(c(0,1,2), c(0,1,2)),
                                  yrange=data.frame(c(0,1,2), c(0,1,2))))  # xrange and yrange have too many columns
   expect_error(sim_thomas_coords(sad2, xmother = xmother, ymother=xmother,
                                  xrange=data.frame(c(0,1), c(0,1), c(0,1)),
                                  yrange=data.frame(c(0,1), c(0,1), c(0,1))))  # xrange and yrange have too many rows
})


test_that("species richness and abundance are correct", {
   sad1 <- sim_sad(100,1000, fix_s_sim = T)
   expect_equal(length(sad1), 100)
   expect_equal(sum(sad1), 1000)

   sim1 <- sim_thomas_community(100, 1000, fix_s_sim = T)
   expect_equal(nrow(sim1$census), 1000)
   expect_equal(length(table(sim1$census$species)), 100)
})

if(FALSE) {
   sad2 <- sim_sad(4, 100)
   mother_points <-  c(1,-1,1,1)
   plot(sim_thomas_coords(sad2, mother_points = mother_points))
   
   xmother <- lapply(1:length(sad2), function(x) runif(2, 0, 1))
   ymother <- lapply(1:length(sad2), function(x) runif(2, 0, 1))
   sim_thomas_coords(sad2, xmother = xmother, ymother= ymother, xrange = c(0, 0.1), yrange= c(0, 0.1))
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
