test_that("Multi Mean Test works", {
  data(iris)
  chart <- iris[, 1:4]
  species <- iris[, 5]

  load("./Expect/meanTestmultiExpected.RData")

  test <- meanTest.multi(chart, species)

  expect_equal(identical(test, output), TRUE)
})
