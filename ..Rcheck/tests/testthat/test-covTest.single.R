test_that("Single Covariance Test works", {
  data(iris)
  X <- iris[, 1:4]

  load("./Expect/covTestsingleExpected.RData")
  test1 <- covTest.single(X, diag(1, 4))
  test2 <- covTest.single(X, diag(1, 4), ball = TRUE)
  test3 <- covTest.single(X, diag(2, 4), ball = TRUE)

  expect_equal(identical(test1, output1), TRUE)
  expect_equal(identical(test2, output2), TRUE)
  expect_equal(identical(test3, output3), TRUE)
})
