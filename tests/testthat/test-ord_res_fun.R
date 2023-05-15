# Phase_space() function
test_that("Phase_space s3.function works", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  dframe <- data.frame(t, conc)
  res <- det_order(dframe)

  expect_s3_class(phase_space(res), "lm")
  expect_equal(phase_space(res), res[[2]])
})

# kin_regr() function
test_that("kin_regr s3.function works", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  dframe <- data.frame(t, conc)
  res <- det_order(dframe)

  expect_s3_class(kin_regr(res), "nls")
  expect_equal(kin_regr(res), res[[4]])
})
# n-order return
test_that("n_value returns the reaction order", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  dframe <- data.frame(t, conc)
  res <- det_order(dframe)

  expect_equal(n_value(res), res[[6]])
  expect_type(n_value(res), "double")
  expect_equal(n_value(res), c(Order = 1))
})

# k_value function
test_that("k_value function returns the k value of nls", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  dframe <- data.frame(t, conc)
  res <- det_order(dframe)

  regr <- nls(y ~ exp(-k * t),
    data = list(
      t = t,
      y = conc
    ),
    start = list(k = 0.17)
  )
  expect_equal(signif(coefficients(regr), 4), signif(k_value(res), 4))
})

test_that("k_value function returns the k value of lm", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- 1 - t * 0.04
  conc <- conc - c(0.01, -0.02, +0.04, -0.01, 0.03, 0.02)
  dframe <- data.frame(t, conc)
  res <- det_order(dframe)

  results(res)
  regr <- lm(conc ~ t)
  expect_equal(signif(abs(coefficients(regr)), 4), signif(k_value(res), 4))
})


test_that("k_value returns nothing if the kinetic regression is
          not significant", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.51, 0.7, 0.12, 0.07)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_message(k_value(res))
})

# results() function
test_that("results works with nls", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  dframe <- data.frame(t, conc)
  res <- det_order(dframe)

  expect_invisible(results(res))
})

test_that("results works with NAN in object 6 (no kinetic order found)", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.51, 0.7, 0.12, 0.07)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_invisible(results(res))
})

test_that("results works with fractional n values", {
  t <- c(0, 4, 8, 12, 16, 20)
  set.seed(1)
  conc <- (0.5 * (0.05 * t + (10^(-0.5)) / 0.5))^-2 + rnorm(length(t), 0, 0.05)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_equal(res[[6]], c(Order = 1.4))
  expect_invisible(results(res))
})

test_that("If estimated model parameters are not significant a particular
          message is displayed - case order==0", {
  t <- c(0, 1, 2, 3, 4, 5)

  set.seed(1)
  conc <- 12 * exp(-0.01 * t) + rnorm(length(t), 0, 1)
  pippo <- data.frame(t, conc, err = abs(rnorm(length(t), 0, 1)))
  expect_no_error(res <- det_order(pippo))
  expect_no_error(results(res))
  expect_equal(results(res), NULL)
})

test_that("If estimated model parameters are not significant a particular
          message is displayed - case order!=0", {
  t <- c(0, 2, 4, 6, 8, 10)
  set.seed(1)
  conc <- 12 * exp(-0.5 * t) + rnorm(length(t), 0, 0.05)
  pippo <- data.frame(t, conc, err = abs(rnorm(length(t), 0, 1)))
  expect_no_error(res <- det_order(pippo))
  res[[3]] <- FALSE
  expect_no_error(results(res))
  expect_equal(results(res), NULL)
})
