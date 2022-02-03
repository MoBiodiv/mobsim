context("Testing patterns of biodiversity functions")

test_that("spec_sample() - gives same rarefied richness as vegan::rarefy", {
   expect_equal(
      as.numeric(spec_sample(abund_vec = 1:100, n = 10)),
      vegan::rarefy(x = 1:100, sample = 10)[1]
   )
})


test_that("rare_curve() - classes and dimensions are correct", {
   S <- 5L
   N <- 50L
   rc1 <- rare_curve(sim_sad(S, N))

   expect_vector(rc1, ptype = numeric(), size = N)
   expect_equal(as.numeric(rc1)[1L], 1L)
   expect_equal(max(rc1), S)
})

test_that("rare_curve() - handles 0 and NA values as expected", {
   sad <- sim_sad(5L, 50L)
   sadNA <- c(NA_integer_, sad)
   sad0 <- c(0L, sad)

   expect_error(rare_curve(sadNA))
   expect_vector(rare_curve(sad0), size = 50L)
})


test_that("div_rect() & div_rand_rect() - return expected value classes and length", {
   S = 5L
   comm1 <- sim_thomas_coords(sim_sad(S, 50L))
   div1 <- div_rect(x0 = 0, y0 = 0, xsize = 1, ysize = 1, comm = comm1)
   divRand <- div_rand_rect(prop_area = 1, n_rect = 2, exclude_zeros = TRUE, comm = comm1)

   expect_vector(div1, ptype = numeric(), size = 6)
   expect_vector(divRand, ptype = numeric(), size = 12)
})

test_that("div_rect() - produces results comparable to vegan", {
   S = 5L
   comm1 <- sim_thomas_coords(sim_sad(S, 50L))
   div1 <- div_rect(x0 = 0, y0 = 0, xsize = 1, ysize = 1, comm = comm1)

   expect_equal(as.numeric(div1[1]), S)
   expect_equal(as.numeric(div1[2]), S)
   expect_equal(as.numeric(div1[3]), vegan::diversity(table(comm1$census$species), "shannon"))
   expect_equal(as.numeric(div1[4]), exp(vegan::diversity(table(comm1$census$species), "shannon")))
   expect_equal(as.numeric(div1[5]), vegan::diversity(table(comm1$census$species), "simpson"))
   expect_equal(as.numeric(div1[6]), 1 / (1 - vegan::diversity(table(comm1$census$species), "simpson")))
})


test_that("div_rand_rect() - produces results comparable to vegan", {
   S = 5L
   comm1 <- sim_thomas_coords(sim_sad(S, 50L))
   divRand <- div_rand_rect(prop_area = 1, n_rect = 2, exclude_zeros = TRUE, comm = comm1)

   expect_equal(as.numeric(divRand[1]), S)
   expect_equal(as.numeric(divRand[3]), S)
   expect_equal(as.numeric(divRand[5]), vegan::diversity(table(comm1$census$species), "shannon"))
   expect_equal(as.numeric(divRand[7]), exp(vegan::diversity(table(comm1$census$species), "shannon")))
   expect_equal(as.numeric(divRand[9]), vegan::diversity(table(comm1$census$species), "simpson"))
   expect_equal(as.numeric(divRand[11]), 1 / (1 - vegan::diversity(table(comm1$census$species), "simpson")))
})


