test_that("single mean test works", {
  data(iris)
  X <- iris[, 1:4]
  mu0 <- c(5.8, 3.0, 4.3, 1.3)

  load("./Expect/meanTestsingleExpected.RData")

  test1 <- meanTest.single(X, mu0)
  test2 <- meanTest.single(X, mu0, Sigma0 = diag(1, 4))

  expect_equal(identical(test1, output1), TRUE)
  expect_equal(identical(test2, output2), TRUE)
})
