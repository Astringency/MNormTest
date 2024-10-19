test_that("Multi Covariance Test works", {
  data(iris)
  chart <- iris[, 1:4]
  species <- iris[, 5]

  load("./Expect/covTestmultiExpected.RData")

  test <- covTest.multi(chart, species)

  expect_equal(identical(test, output), TRUE)
})
