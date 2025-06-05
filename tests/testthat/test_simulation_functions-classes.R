test_that("classes are correct", {
  sad1 <- sim_sad(s_pool = 3L, n_sim = 10L)
  expect_s3_class(sad1, "sad")
  expect_s3_class(sad1, "integer")

  xmother <- lapply(1:length(sad1), function(x) runif(1, 0, 1))
  ymother <- lapply(1:length(sad1), function(x) runif(1, 0, 1))

  expect_s3_class(sim_poisson_coords(sad1), "community")
  expect_s3_class(sim_poisson_community(3L, 10L), "community")
  expect_s3_class(sim_thomas_coords(sad1), "community")
  expect_s3_class(sim_thomas_coords(sad1, mother_points = 0), "community")
  expect_s3_class(sim_thomas_coords(sad1, mother_points = 100), "community")
  expect_s3_class(sim_thomas_coords(sad1, mother_points = 4.2), "community")
  expect_s3_class(
    sim_thomas_coords(sad1, xmother = xmother, ymother = ymother),
    "community"
  )
  expect_s3_class(sim_thomas_coords(sad1), "community")
  expect_s3_class(sim_thomas_community(3L, 10L), "community")
  expect_s3_class(
    community(runif(10L), runif(10L), rep("specA", 10L)),
    "community"
  )
  expect_s3_class(community_to_sad(sim_poisson_community(3L, 10L)), "sad")
})
