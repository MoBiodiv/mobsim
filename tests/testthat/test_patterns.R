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
