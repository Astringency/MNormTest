pkgname <- "MNormTest"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('MNormTest')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("covTest.multi")
### * covTest.multi

flush(stderr()); flush(stdout())

### Name: covTest.multi
### Title: Multiple Covariance Matrix Hypothesis Testing
### Aliases: covTest.multi

### ** Examples

data(iris)
chart <- iris[, 1:4]
species <- iris[, 5]
# carry out the test
test1 <- covTest.multi(chart, species)
test2 <- covTest.multi(chart, species, detail = FALSE)
# get the elements
test1$Stat
test1$SampMeanT
test1$sampleSize



cleanEx()
nameEx("covTest.single")
### * covTest.single

flush(stderr()); flush(stdout())

### Name: covTest.single
### Title: Single Covariance Matrix Hypothesis Testing
### Aliases: covTest.single

### ** Examples

data(iris)
X <- iris[, 1:4]
# carry out the test
test1 <- covTest.single(X, diag(1, 4))
test2 <- covTest.single(X, diag(1, 4), ball = TRUE)
test3 <- covTest.single(X, diag(2, 4), ball = TRUE)
test4 <- covTest.single(X, diag(1, 4), detail = FALSE)
# get the elements
test1$Stat
test2$Df
test3$sigma.hat



cleanEx()
nameEx("indTest.multi")
### * indTest.multi

flush(stderr()); flush(stdout())

### Name: indTest.multi
### Title: Multivariate Normal Independence Test
### Aliases: indTest.multi

### ** Examples

data(iris)
chart <- iris[, 1:4]
# carry out the test
test1 <- indTest.multi(chart)
test2 <- indTest.multi(chart, subdim = c(2, 1, 1))
test3 <- indTest.multi(chart, detail = FALSE)
# get the elements
test1$Stat
test1$SampMean
test2$SampAii



cleanEx()
nameEx("meanTest.multi")
### * meanTest.multi

flush(stderr()); flush(stdout())

### Name: meanTest.multi
### Title: Multiple Mean Vectors Hypothesis Testing
### Aliases: meanTest.multi

### ** Examples

data(iris)
chart <- iris[, 1:4]
species <- iris[, 5]
# carry out the test
test1 <- meanTest.multi(chart, species)
test2 <- meanTest.multi(chart, species, detail = FALSE)
# get the elements
test1$Stat
test1$SampMeanT
test1$sampleSize



cleanEx()
nameEx("meanTest.single")
### * meanTest.single

flush(stderr()); flush(stdout())

### Name: meanTest.single
### Title: Single Mean Vector Hypothesis Testing
### Aliases: meanTest.single

### ** Examples

data(iris)
X <- iris[, 1:4]
mu0 <- c(5.8, 3.0, 4.3, 1.3)
# carry out the test
test1 <- meanTest.single(X, mu0)
test2 <- meanTest.single(X, mu0, Sigma0 = diag(1, 4))
test3 <- meanTest.single(X, mu0, detail = FALSE)
# get the elements
test1$Stat
test1$SampMean
test1$SampA
test1$Df



cleanEx()
nameEx("meanTest.two")
### * meanTest.two

flush(stderr()); flush(stdout())

### Name: meanTest.two
### Title: Two Mean Vectors Hypothesis Testing
### Aliases: meanTest.two

### ** Examples

data(iris)
X <- iris[1:50, 1:4]
Y <- iris[51:100, 1:4]
# carry out the test
test1 <- meanTest.two(X, Y)
test2 <- meanTest.two(X, Y, detail = TRUE)
test3 <- meanTest.two(X, Y, equal = FALSE, method = "Coupled")
test4 <- meanTest.two(X, Y, equal = FALSE, method = "Transformed")
# get the elements
test1$Stat
test1$SampMean1
test3$SampMeanC
test4$dataT



cleanEx()
nameEx("meancov.Test")
### * meancov.Test

flush(stderr()); flush(stdout())

### Name: meancov.Test
### Title: Mean and Covariance Matrix Hypothesis Testing (Simultaneously)
### Aliases: meancov.Test

### ** Examples

data(iris)
chart <- iris[, 1:4]
species <- iris[, 5]
# carry out the test
test1 <- meancov.Test(chart, species)
test2 <- meancov.Test(chart, species, detail = FALSE)
# get the elements
test1$Stat
test1$SampMeanT



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
