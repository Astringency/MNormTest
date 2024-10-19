test_that("Mean Covariance Test works", {
  data(iris)
  chart <- iris[, 1:4]
  species <- iris[, 5]

  load("./Expect/meancovTestExpected.RData")

  test <- meancov.Test(chart, species)

  expect_equal(identical(test, output), TRUE)
})
