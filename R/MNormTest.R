#' @title Test of Single Mean Vector
#' @description Hypothesis testing of the mean vector of a single multivariate normal population when the covariance matrix is known or unknown.
#' @author Xifeng Zhang
#' @param data The data matrix which is a matrix or data frame.
#' @param mu0 The mean vector when the null hypothesis is true.
#' @param Sigma0 (optional) The population covariance matrix. Default is FALSE which means the covariance matrix is unknown.
#' @param alpha The significance level. Default is 0.05.
#' @import stats
#' @import utils
#' @references 高惠璇. 应用多元统计分析. 北京大学出版社, 2005: 66-68.
#' @return A list containing the hypothesis, sample mean, sample diviation, statistics, df of T2, df of F, p value, critical value and conclusion.
#' @export
mean.test.single <- function(data, mu0, Sigma0 = FALSE, alpha = 0.05) {
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)

  X.bar <- round(apply(data, 2, mean), 4)

  if (isFALSE(Sigma0)) {
    A <- (n - 1) * cov(data)
    T2 <- (n - 1) * n * t(X.bar - mu0) %*% solve(A) %*% (X.bar - mu0)
    F <- (n - p) / ((n - 1) * p) * T2
    p.value <- min(pf(F, p, n - p), 1 - pf(F, p, n - p))
    critical.value <- qf(1 - alpha, p, n - p)
    reject <- ifelse(F > critical.value, "Reject", "Not Reject")
    statistics <- data.frame("Statistics" = c("Hotelling T2", "F"), "Value" = c(T2, F))
    return(list(
      "Hypothesis" = paste("H0: mu = (", toString(mu0), ") when Sigma is unknown"),
      "Sample Mean" = X.bar,
      "Sample Diviation" = A,
      "Statistics" = statistics,
      "df of T2" = n - 1,
      "df of F" = paste("df1 = ", p, ", df2 = ", n - p),
      "p value" = p.value,
      "Critical Value" = critical.value,
      "conclusion" = reject
    ))
  } else {
    T0 <- n * t(X.bar - mu0) %*% solve(Sigma0) %*% (X.bar - mu0)
    p.value <- min(pchisq(T0, p), 1 - pchisq(T0, p))
    critical.value <- qchisq(1 - alpha, p)
    reject <- ifelse(T0 > critical.value, "Reject", "Not Reject")
    return(list(
      "Hypothesis" = paste("H0: mu = (", toString(mu0), ") when Sigma0 is known"),
      "Sample Mean" = X.bar,
      "Population Covariance" = Sigma0,
      "Statistics" = c("T0" = T0),
      "df of T0" = p,
      "p value" = p.value,
      "Critical Value" = critical.value,
      "conclusion" = reject
    ))
  }
}

#' @title Test of Two Mean Vectors
#' @description Hypothesis testing for equality of the mean vectors of two multivariate normal populations when the covariance matrices are equal or unequal.
#' @param data1 A matrix or data frame of group 1.
#' @param data2 A matrix or data frame of group 2.
#' @param alpha The significance level. Default is 0.05.
#' @param equal A boolean value. Default is TRUE. If TRUE, the covariance matrix is equal. If FALSE, the covariance matrix is not equal.
#' @param method A string value. Default is "None". When equal is FALSE, you must choose a method in "Coupled" or "Transformed". Choose "Coupled" when the sample size of two groups is equal. Choose "Transformed" when the sample size of two groups is not equal. If you want to use the likelihood ratio test, please use the function "mean.test.multi".
#' @import stats
#' @import utils
#' @references 高惠璇. 应用多元统计分析. 北京大学出版社, 2005: 76-80.
#' @return A list containing the hypothesis, sample mean, sample diviation, statistics, df of T2, df of F, p value, critical value and conclusion.
#' @export
mean.test.two <- function(data1, data2, alpha = 0.05, equal = TRUE, method = c("None", "Coupled", "Transformed")) {
  data1 <- as.matrix(data1)
  data2 <- as.matrix(data2)

  n1 <- nrow(data1)
  n2 <- nrow(data2)
  p <- ncol(data1)

  if (equal) {
    X.bar <- apply(data1, 2, mean)
    A1 <- (n1 - 1) * cov(data1)
    Y.bar <- apply(data2, 2, mean)
    A2 <- (n2 - 1) * cov(data2)
    A <- (A1 + A2) / (n1 + n2 - 2)

    T2 <- (n1 * n2 / (n1 + n2)) * t(X.bar - Y.bar) %*% solve(A) %*% (X.bar - Y.bar)
    F <- (n1 + n2 - p - 1) / ((n1 + n2 - 2) * p) * T2

    Critical.Value <- qf(1 - alpha, p, n1 + n2 - p - 1)
    Reject <- ifelse(F > Critical.Value, "Reject", "Not Reject")
    p.value <- min(pf(F, p, n1 + n2 - p - 1), 1 - pf(F, p, n1 + n2 - p - 1))

    return(list(
      "Hypothesis" = "H0: mu1 = mu2, with unknown but equal covariance matrix",
      "Sample Mean X" = X.bar,
      "Sample Mean Y" = Y.bar,
      "Sample Diviation A1" = A1,
      "Sample Diviation A2" = A2,
      "Statistics" = c("Hotelling T2" = T2, "F" = F),
      "df of T2" = n1 + n2 - 2,
      "df of F" = c("df1" = p, "df2" = n1 + n2 - p - 1),
      "p value" = p.value,
      "Critical Value" = Critical.Value,
      "conclusion" = Reject
    ))
  } else {
    if (method == "None") {
      print("Please choose a method in 'Coupled'、'Transformed' or 'LRT'")
    } else if (method == "Coupled" && n1 == n2) {
      datacoupled <- data1 - data2
      X.bar <- apply(datacoupled, 2, mean)
      A <- (n1 - 1) * cov(datacoupled)
      T2 <- (n1 - 1) * t(X.bar) %*% solve(A) %*% X.bar
      F <- (n1 - p) / ((n1 - 1) * p) * T2
      Critical.Value <- qf(1 - alpha, p, n1 - p)
      Reject <- ifelse(F > Critical.Value, "Reject", "Not Reject")
      p.value <- min(pf(F, p, n1 - p), 1 - pf(F, p, n1 - p))
      return(list(
        "Hypothesis" = paste("H0: mu1 = mu2, with unknown but different covariance matrix"),
        "Method" = "Construct a new sample by subtracting two samples",
        "Sample Mean" = X.bar,
        "Sample Diviation" = A,
        "Statistics" = c("Hotelling T2" = T2, "F" = F),
        "df of T2" = n1 - 1,
        "df of F" = c("df1" = p, "df2" = n1 - p),
        "p value" = p.value,
        "Critical Value" = Critical.Value,
        "conclusion" = Reject
      ))
    } else if (method == "Coupled" && n1 != n2) {
      print("The sample size of two groups should be equal when using Coupled method!")
    } else if (method == "Transformed") {
      data.z <- data1 - sqrt(n1 / n2) * data2 + 1 / sqrt(n1 * n2) * apply(data2[1:n1, ], 2, sum) - 1 / n2 * apply(data2, 2, sum)
      Z.bar <- apply(data.z, 2, mean)
      A <- (n1 - 1) * cov(data.z)
      T2 <- (n1 - 1) * t(Z.bar) %*% solve(A) %*% Z.bar
      F <- (n1 - p) / ((n1 - 1) * p) * T2
      Critical.Value <- qf(1 - alpha, p, n1 - p)
      Reject <- ifelse(F > Critical.Value, "Reject", "Not Reject")
      p.value <- min(pf(F, p, n1 - p), 1 - pf(F, p, n1 - p))
      return(list(
        "Hypothesis" = paste("H0: mu1 = mu2, with unknown but different covariance matrix"),
        "Method" = "Transformed method",
        "Sample Mean" = Z.bar,
        "Sample Diviation" = A,
        "Statistics" = c("T2" = T2, "F" = F),
        "df of T2" = n1 - 1,
        "df of F" = c("df1" = p, "df2" = n1 - p),
        "p value" = p.value,
        "Critical Value" = Critical.Value,
        "conclusion" = Reject
      ))
    }
  }
}

#' @title Test of Multiple Mean Vectors
#' @description Mean vector test for multiple multivariate normal totals when the covariance array is equal (multivariate analysis of variance).
#' @param X The data matrix which is a matrix or data frame.
#' @param label A vector of group labels.
#' @param alpha The significance level. Default is 0.05.
#' @import stats
#' @import utils
#' @references 高惠璇. 应用多元统计分析. 北京大学出版社, 2005: 80-83.
#' @return A list containing the hypothesis, sample size, total sample mean, within group mean, total sum of squares, within sum of squares, between sum of squares, statistics, df of Wilk's Lambda, df of Bartlett's Chi2, df of Rao's F, p value, approximate conclusion.
#' @export
mean.test.multi <- function(X, label, alpha = 0.05) {
  data <- cbind(X, label)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  k <- length(unique(label))

  X.bar <- apply(X, 2, mean)
  T <- (n - 1) * cov(X)
  data.split <- split(data[, 1:p], data[, p + 1])
  Within.mean <- lapply(data.split, function(x) {
    apply(x, 2, mean)
  })
  At <- lapply(data.split, function(x) {
    (nrow(x) - 1) * cov(x)
  })
  A <- Reduce("+", At)
  nt <- sapply(data.split, nrow)

  Lambda <- det(A) / det(T)
  df.Lambda <- c("df1" = p, "df2" = n - k, "df3" = k - 1)

  Bartlett.chi2 <- -(df.Lambda[2] + df.Lambda[3] - (df.Lambda[1] + df.Lambda[3] + 1) / 2) * log(Lambda)
  df.chi2 <- df.Lambda[1] * df.Lambda[3]
  p.value.chi2 <- min(pchisq(Bartlett.chi2, df.chi2), 1 - pchisq(Bartlett.chi2, df.chi2))

  t <- df.Lambda[2] + df.Lambda[3] - (df.Lambda[1] + df.Lambda[3] + 1) / 2
  s <- sqrt((df.Lambda[1]^2 * df.Lambda[3]^2 - 4) / (df.Lambda[1]^2 + df.Lambda[3]^2 - 5))
  lambda <- (df.Lambda[1] * df.Lambda[3] - 2) / 4
  Rao.F <- (1 - Lambda^(1 / s)) * (t * s - 2 * lambda) / (Lambda^(1 / s)) * (df.Lambda[1] * df.Lambda[3])
  df.F <- c("df1" = df.Lambda[1] * df.Lambda[3], "df2" = round(t * s - 2 * lambda))
  p.value.F <- min(pf(Rao.F, df.F[1], df.F[2]), 1 - pf(Rao.F, df.F[1], df.F[2]))

  Reject <- c(
    "Bartlett" = ifelse(p.value.chi2 < alpha, "Reject", "Not Reject"),
    "Rao" = ifelse(p.value.F < alpha, "Reject", "Not Reject")
  )

  return(list(
    "Hypothesis" = paste("H0: mu1 = mu2 = ... = muk, k = ", k),
    "Sample Size" = nt,
    "Total Sample Mean" = X.bar,
    "Within Group Mean" = Within.mean,
    "Total Sum of Squares" = T,
    "Within Sum of Squares (Total)" = A,
    "Within Sum of Squares" = At,
    "Between Sum of Squares" = T - A,
    "Statistics" = c("Wilk's Lambda" = Lambda, "Bartlett's Chi2" = Bartlett.chi2, "Rao's F" = Rao.F),
    "df of Wilk's Lambda" = df.Lambda,
    "df of Bartlett's Chi2" = df.chi2,
    "df of Rao's F" = df.F,
    "p value" = c("Bartlett's Chi2" = p.value.chi2, "Rao's F" = p.value.F),
    "Approx Reject" = Reject
  ))
}

#' @title Test of Single Covariance Matrix
#' @description Hypothesis testing of covariance matrices for a single normal population when the mean vector is unknown.
#' @param data The data matrix which is a matrix or data frame.
#' @param Sigma0 The covariance matrix when the null hypothesis is true.
#' @param ball A boolean value. Default is FALSE. If FALSE, the covariance matrix is Sigma0 (known). If TRUE and the Sigma0 unit matrix, the Mauchly ball test is performed. If TRUE but Sigma0 is not a unit matrix, the covariance array is tested to see if it is sigma^2*Sigma0 (sigma^2 is unknown).
#' @param alpha The significance level. Default is 0.05.
#' @import stats
#' @import utils
#' @import Rmpfr
#' @references 高惠璇. 应用多元统计分析. 北京大学出版社, 2005: 83-88.
#' @return A list containing the hypothesis, sample mean, sample diviation, statistics, df of Chi2, p value, critical value and conclusion.
#' @export
cov.test.single <- function(data, Sigma0, ball = FALSE, alpha = 0.05) {
  if (require("Rmpfr")) {
    print("Successfully loaded package 'Rmpfr'.")
  } else {
    print("Lack of package 'Rmpfr'. Now installing...")
    install.packages("Rmpfr")
    if (require("Rmpfr")) {
      print("Successfully installed package 'Rmpfr'.")
    } else {
      stop("Failed to install package 'Rmpfr'.")
    }
  }

  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)

  X.bar <- round(apply(data, 2, mean), 4)
  A <- (n - 1) * cov(data)

  if (isFALSE(ball)) {
    lambda <- exp(-1 / 2 * trace(A %*% solve(Sigma0))) * mpfr(det(A %*% solve(Sigma0)), 256)^(n / 2) * exp(n * p / 2) * n^(-np / 2)
    chi2 <- -2 * log(lambda)

    lambda <- as.numeric(lambda)
    chi2 <- as.numeric(chi2)

    df.chi2 <- p * (p + 1) / 2
    p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
    critical.value <- qchisq(1 - alpha, df.chi2)
    reject <- ifelse(p.value < alpha, "Reject", "Not Reject")
    statistics <- data.frame("Statistics" = c("Likelihood Ratio", "Chi2"), "Value" = c(lambda, chi2))

    return(list(
      "Hypothesis" = "H0: Sigma = Sigma0",
      "Sample Mean" = X.bar,
      "Sample Diviation" = A,
      "Statistics" = statistics,
      "df of Chi2" = df.chi2,
      "p value" = p.value,
      "Critical Value" = critical.value,
      "conclusion" = reject
    ))
  } else {
    sigma.hat <- trace(solve(Sigma0) %*% A) / (n * p)
    lambda <- (mpfr(det(solve(Sigma0) %*% A), 256)^(n / 2)) / ((trace(solve(Sigma0) %*% A) / p)^(n * p / 2))
    W <- lambda^(2 / n)
    chi2 <- -((n - 1) - (2 * p^2 + p + 2) / (6 * p)) * log(W)

    lambda <- as.numeric(lambda)
    W <- as.numeric(W)
    chi2 <- as.numeric(chi2)

    df.chi2 <- p * (p + 1) / 2 - 1
    p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
    critical.value <- qchisq(1 - alpha, df.chi2)
    reject <- ifelse(p.value < alpha, "Reject", "Not Reject")
    statistics <- data.frame("Statistics" = c("Likelihood Ratio", "W", "Chi2"), "Value" = c(lambda, W, chi2))

    return(list(
      "Hypothesis" = ifelse((Sigma0 == diag(1, p)), "Mauchly’s test of sphericity. H0: Sigma = sigma2 * Ip", "H0: Sigma = sigma2 * Sigma0"),
      "Sample Mean" = X.bar,
      "Sample Diviation" = A,
      "sigma.hat" = sigma.hat,
      "Statistics" = statistics,
      "df of Chi2" = df.chi2,
      "p value" = p.value,
      "Critical Value" = critical.value,
      "conclusion" = reject
    ))
  }
}

#' @title Test of Multiple Covariance Matrix
#' @description Tests of multiple multivariate normal overall covariance arrays when the overall mean is unknown.
#' @param X The data matrix which is a matrix or data frame.
#' @param label A vector of group labels.
#' @param alpha The significance level. Default is 0.05.
#' @import stats
#' @import utils
#' @import Rmpfr
#' @references 高惠璇. 应用多元统计分析. 北京大学出版社, 2005: 88-89.
#' @return A list containing the hypothesis, sample size, total sample mean, within group mean, within group sample covariance, total sum of squares, within sum of squares, statistics, modify factor, df of Chi2, p value, critical value and conclusion.
#' @export
cov.test.multi <- function(X, label, alpha = 0.05) {
  if (require("Rmpfr")) {
    print("Successfully loaded package 'Rmpfr'.")
  } else {
    print("Lack of package 'Rmpfr'. Now installing...")
    install.packages("Rmpfr")
    if (require("Rmpfr")) {
      print("Successfully installed package 'Rmpfr'.")
    } else {
      stop("Failed to install package 'Rmpfr'.")
    }
  }

  data <- data.frame(cbind(X, label))
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  k <- length(unique(label))

  X.bar <- apply(X, 2, mean)

  data.split <- split(data[, 1:p], data[, p + 1])
  nt <- sapply(data.split, nrow)
  Within.mean <- lapply(data.split, function(x) {
    apply(x, 2, mean)
  })
  At <- lapply(data.split, function(x) {
    (nrow(x) - 1) * cov(x)
  })
  A <- Reduce("+", At)

  St <- lapply(data.split, function(x) {
    cov(x)
  })
  St.root <- lapply(data.split, function(x) {
    mpfr(det(cov(x)), 256)^(-(nrow(x) - 1) / 2)
  })
  St.sum <- lapply(data.split, function(x) {
    (nrow(x) - 1) * log(mpfr(det(cov(x)), 256))
  })
  St.star <- lapply(data.split, function(x) {
    cov(x) * (nrow(x) - 1) / nrow(x)
  })
  St.star.root <- lapply(data.split, function(x) {
    mpfr(det(cov(x) * (nrow(x) - 1) / nrow(x)), 256)^(-nrow(x) / 2)
  })

  lambda <- det(A / nrow(data))^(-n / 2) / Reduce("*", St.star.root)
  lambda.star <- det(A / (n - k))^(-(n - k) / 2) / Reduce("*", St.root)
  M <- -2 * log(lambda.star)
  M.star <- (n - k) * log(det(A / (n - k))) - Reduce("+", St.sum)

  lambda <- as.numeric(lambda)
  lambda.star <- as.numeric(lambda.star)
  M <- as.numeric(M)
  M.star <- as.numeric(M.star)

  if (identical(nt, rep(1, k)) == TRUE) {
    d <- (2 * p^2 + 3 * p - 1) * (k + 1) / (6 * (p + 1) * (n - k))
  } else {
    d <- (2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (k - 1)) * (sum(1 / (nt - 1)) - 1 / (n - k))
  }

  chi2 <- (1 - d) * M
  df.chi2 <- p * (p + 1) * (k - 1) / 2
  p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
  critical.value <- qchisq(1 - alpha, df.chi2)
  reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

  statistics <- data.frame(
    "Statistics" = c("Likelihood Ratio", "Modified Likelihood Ratio", "M", "Chi2"),
    "Value" = c(lambda, lambda.star, M.star, chi2)
  )

  return(list(
    "Hypothesis" = paste("H0: Sigma1 = Sigma2 = ... = Sigmak, k = ", k),
    "Sample Size" = nt,
    "Total Sample Mean" = X.bar,
    "Within Group Mean" = Within.mean,
    "Within Group Sample Covariance" = St,
    "Total Sum of Squares" = A,
    "Within Sum of Squares" = At,
    "Statistics" = statistics,
    "Modify Factor" = d,
    "df of Chi2" = df.chi2,
    "p value" = p.value,
    "Critical Value" = critical.value,
    "conclusion" = reject
  ))
}

#' @title Test of Mean and Covariance Matrix at the Same Time
#' @description Simultaneous testing of multiple multivariate normal overall mean vectors and covariance matrices.
#' @param X The data matrix which is a matrix or data frame.
#' @param label A vector of group labels.
#' @param alpha The significance level. Default is 0.05.
#' @import stats
#' @import utils
#' @import Rmpfr
#' @references 高惠璇. 应用多元统计分析. 北京大学出版社, 2005: 90-91.
#' @return A list containing the hypothesis, sample size, total sample mean, within group mean, total sum of squares, within sum of squares, statistics, modify factor, df of Chi2, p value, critical value and conclusion.
#' @export
MeanCov.test <- function(X, label, alpha = 0.05) {
  if (require("Rmpfr")) {
    print("Successfully loaded package 'Rmpfr'.")
  } else {
    print("Lack of package 'Rmpfr'. Now installing...")
    install.packages("Rmpfr")
    if (require("Rmpfr")) {
      print("Successfully installed package 'Rmpfr'.")
    } else {
      stop("Failed to install package 'Rmpfr'.")
    }
  }

  data <- data.frame(cbind(X, label))
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  k <- length(unique(label))

  X.bar <- apply(X, 2, mean)
  data.split <- split(data[, 1:p], data[, p + 1])
  nt <- sapply(data.split, nrow)
  Within.mean <- lapply(data.split, function(x) {
    apply(x, 2, mean)
  })
  At <- lapply(data.split, function(x) {
    (nrow(x) - 1) * cov(x)
  })
  A <- Reduce("+", At)
  T <- (n - 1) * cov(X)

  lambda <- Reduce("*", lapply(data.split, function(x) {
    mpfr(det((nrow(x) - 1) * cov(x)), 256)^(nrow(x) / 2)
  })) / mpfr(det(T), 256)^(n / 2) * n^(n * p / 2) / mpfr(Reduce("*", lapply(data.split, function(x) {
    nrow(x)^(nrow(x) * p / 2)
  })), 256)

  lambda.star <- Reduce("*", lapply(data.split, function(x) {
    mpfr(det((nrow(x) - 1) * cov(x)), 256)^((nrow(x) - 1) / 2)
  })) / mpfr(det(T), 256)^((n - k) / 2) * (n - k)^((n - k) * p / 2) / mpfr(Reduce("*", lapply(data.split, function(x) {
    (nrow(x) - 1)^((nrow(x) - 1) * p / 2)
  })), 256)

  d <- (Reduce("+", lapply(data.split, function(x) {
    sum(1 / (nrow(x) - 1))
  })) - sum(1 / (n - k))) * ((2 * p^2 + 3 * p - 1) / (6 * (p + 3) * (k - 1))) - (p - k + 2) / ((n - k) * (p + 3))

  M <- -2 * log(lambda.star)
  chi2 <- (1 - d) * M

  df.chi2 <- p * (p + 3) * (k - 1) / 2

  lambda <- as.numeric(lambda)
  lambda.star <- as.numeric(lambda.star)
  M <- as.numeric(M)
  chi2 <- as.numeric(chi2)

  p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
  critical.value <- qchisq(1 - alpha, df.chi2)
  reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

  statistics <- data.frame("Statistics" = c("Likelihood Ratio", "Likelihood Ratio (Modified)", "M", "Chi2"), "Value" = c(lambda, lambda.star, M, chi2))

  return(list(
    "Hypothesis" = paste("H0: mu1 = mu2 = ... = muk, Sigma1 = Sigma2 = ... = Sigmak, k = ", k),
    "Sample Size" = nt,
    "Total Sample Mean" = X.bar,
    "Within Group Mean" = Within.mean,
    "Total Sum of Squares" = T,
    "Within Sum of Squares" = A,
    "Statistics" = statistics,
    "Modify Factor" = d,
    "df of Chi2" = df.chi2,
    "p value" = p.value,
    "Critical Value" = critical.value,
    "conclusion" = reject
  ))
}

#' @title Test of Independence
#' @description Independence tests for multivariate normal populations.
#' @param data The data matrix which is a matrix or data frame.
#' @param subdim The dimensions of submatrices. The default is FALSE, which means the independence of all components of the random vector will be tested.
#' @param alpha The significance level. Default is 0.05.
#' @import stats
#' @import utils
#' @references 高惠璇. 应用多元统计分析. 北京大学出版社, 2005: 92-94.
#' @return A list containing the hypothesis, dimension, sample mean, sample diviation, sample diviation of submatrix, statistics, modify factor, df of Chi2, p value, critical value and conclusion.
#' @export
ind.test.multi <- function(data, subdim = FALSE, alpha = 0.05) {
  n <- nrow(data)
  p <- ncol(data)

  X.bar <- round(apply(data, 2, mean), 4)
  A <- (n - 1) * cov(data)
  Aii <- list()
  subdata <- list()

  if (isFALSE(subdim)) {
    subdim <- rep(1, p)
  } else {
    subdim <- subdim
  }

  subcol <- c(0, cumsum(subdim))

  for (i in 1:length(subdim)) {
    datai <- data[, (subcol[i] + 1):subcol[i + 1]]
    Aii <- append(Aii, list((n - 1) * cov(as.matrix(datai))))
  }

  V <- (det(A) / Reduce("*", lapply(Aii, det)))
  lambda <- V^(n / 2)

  b <- n - 3 / 2 - (p^3 - sum(subdim^3)) / (3 * (p^2 - sum(subdim^2)))
  chi2 <- -b * log(V)

  df.chi2 <- 1 / 2 * (p * (p + 1) - sum(subdim * (subdim + 1)))

  p.value <- min(pchisq(chi2, df.chi2), 1 - pchisq(chi2, df.chi2))
  critical.value <- qchisq(1 - alpha, df.chi2)
  reject <- ifelse(p.value < alpha, "Reject", "Not Reject")

  statistics <- data.frame("Statistics" = c("V", "Likelihood Ratio", "Chi2"), "Value" = c(V, lambda, chi2))

  return(list(
    "Hypothesis" = "H0: The components are independent",
    "Dimension" = subdim,
    "Sample Mean" = X.bar,
    "Sample Diviation" = A,
    "Sample Diviation of Submatrix" = Aii,
    "Statistics" = statistics,
    "Modify Factor" = b,
    "df of Chi2" = df.chi2,
    "p value" = p.value,
    "Critical Value" = critical.value,
    "conclusion" = reject
  ))
}
