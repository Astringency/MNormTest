test_that("independent Test works", {
  data(iris)
  chart <- iris[, 1:4]

  load("./Expect/indTestmultiExpected.RData")

  test1 <- indTest.multi(chart)
  test2 <- indTest.multi(chart, subdim = c(2, 1, 1))

  expect_equal(identical(output1, test1), TRUE)
  expect_equal(identical(output2, test2), TRUE)
})
