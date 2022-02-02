context("Testing patterns of biodiversity functions")

test_that("spec_sample gives same rarefied richness as vegan::rarefy", {
   expect_equal(
      as.numeric(spec_sample(abund_vec = 1:100, n = 10)),
      vegan::rarefy(x = 1:100, sample = 10)[1]
   )
})
