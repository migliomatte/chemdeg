test_that("plot_ord works with weighted data", {
  t <- c(0, 4, 8, 12, 16, 20)
  set.seed(1)
  conc <- (0.5 * (0.05 * t + (10^(-0.5)) / 0.5))^-2 + rnorm(length(t), 0, 0.05)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_no_error(plot_ord(res))
})

test_that("plot_ord works with unweighted data", {
  t <- c(0, 4, 8, 12, 16, 20)
  set.seed(1)
  conc <- (0.5 * (0.05 * t + (10^(-0.5)) / 0.5))^-2 + rnorm(length(t), 0, 0.05)
  dframe <- data.frame(t, conc)
  res <- det_order(dframe)
  expect_no_error(plot_ord(res))
})

test_that("plot_ord does not give error with no reaction order determined", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.51, 0.7, 0.12, 0.07)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_equal(res[[6]], NaN, ignore_attr = TRUE)
  expect_no_error(plot_ord(res))
})

test_that("plot_ord works when 0-order model is determined", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 0.80, 0.69, 0.59, 0.51)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_no_error(plot_ord(res))
})

test_that("plot_ord works only with ord_res objects", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.51, 0.7, 0.12, 0.07)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  pippo <- lm(conc ~ t)
  pippo2 <- nls(conc ~ exp(-k * t), start = list(k = 0.5), weights = 1 / err^2)
  dframe <- data.frame(t, conc, err)
  expect_error(plot_ord(t))
  expect_error(plot_ord(dframe))
  expect_error(plot_ord(pippo))
  expect_error(plot_ord(pippo2))
})

test_that("If estimated model parameters are not significant a particular
          message is displayed - case order==0", {
  t <- c(0, 1, 2, 3, 4, 5)

  set.seed(1)
  conc <- 12 * exp(-0.01 * t) + rnorm(length(t), 0, 1)
  pippo <- data.frame(t, conc, err = abs(rnorm(length(t), 0, 1)))
  expect_no_error(res <- det_order(pippo))
  expect_no_error(plot_ord(res))
  expect_message(
    plot_ord(res),
    "NB: parameters of the regression are not significant!"
  )
})
