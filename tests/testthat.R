library(testthat)
library(mobsim)

#test_check("mobsim")
#test_package("mobsim", reporter="summary")
test_dir("./tests/testthat", reporter="summary")
