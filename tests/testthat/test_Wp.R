test_that("Wp function works",
  expect_type(Wp(c(1,2,3), c(1/3, 1/3, 1/3), c(2,3,4), c(1/4, 1/4, 1/2), 1), "double")
)
