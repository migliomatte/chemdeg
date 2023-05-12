# lmp
test_that("lmp works with lm weighted objects", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 0.80, 0.69, 0.59, 0.51)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  pippo <- lm(conc ~ t, weights = 1 / err^2)

  expect_equal(lmp(pippo), 5.425e-07, tolerance = 1e-05)
})

test_that("lmp works with lm unweighted objects", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 0.80, 0.69, 0.59, 0.51)
  pippo <- lm(conc ~ t)
  summary(pippo)
  expect_equal(lmp(pippo), 1.854e-06, tolerance = 1e-05)
})

test_that("lmp returns error with non lm obkects", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 0.80, 0.69, 0.59, 0.51)
  dframe <- data.frame(t, conc)
  pippo <- nls(conc ~ k * t^2, start = list(k = 0.5))
  expect_error(lmp(pippo))
  expect_error(lmp(t))
  expect_error(lmp(dframe))
})

# f_gen
test_that("f_gen works with n=1", {
  expect_type(f_gen(1), "language")
  expect_equal(f_gen(1), as.formula("y ~ y0 * exp(-k * t)"), ignore_attr = TRUE)
})

test_that("f_gen works with n!=1", {
  expect_type(f_gen(7), "language")
  expect_equal(f_gen(2), as.formula("y ~ (1 * (k * t + (y0^(-1))/1))^-1"),
    ignore_attr = TRUE
  )
})

# dup_rem
test_that("dup_rem removes duplicates in dataframe with 2 columns", {
  pippo <- data.frame(t = c(0, 1, 2, 3, 4, 5, 6, 7),
                      x = c(10, 10, 8, 7, 7, 7, 6, 6))
  res <- data.frame(t = c(0.5, 2, 4.25, 6.5), x = c(10, 8, 7, 6))
  expect_equal(dup_rem(pippo), res, ignore_attr = TRUE)
})

test_that("dup_rem removes duplicates in dataframe with 3 columns", {
  pippo <- data.frame(
    t = c(0, 1, 2, 3, 4, 5, 6, 7),
    x = c(10, 10, 8, 7, 7, 7, 6, 6),
    err = c(0.5, 0.6, 0.4, 0.7, 0.8, 0.6, 0.12, 0.24)
  )
  res <- data.frame(t = c(0.5, 2, 4.25, 6.5),
                    x = c(10, 8, 7, 6),
                    err = c(1.1, 0.4, 2.1, 0.36))
  expect_equal(dup_rem(pippo), res, ignore_attr = TRUE)
})

# check_fit
test_that("check fit returns true with a significant linear regression", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 0.80, 0.69, 0.59, 0.51)
  pippo <- lm(conc ~ t)

  expect_lt(summary(pippo)$coefficients[2, 4], 0.05)
  expect_equal(check_fit(pippo), TRUE, ignore_attr = TRUE)
})

test_that("check fit returns false with a non-significant linear regression", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 1.2, 1.1, 1, 0.96)
  pippo <- lm(conc ~ t)

  expect_gt(summary(pippo)$coefficients[2, 4], 0.05)
  expect_equal(check_fit(pippo), FALSE, ignore_attr = TRUE)
})

test_that("check fit returns true with a significant linear regression", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.26, 0.12, 0.6, 0.3)
  pippo <- nls(conc ~ exp(-k * t), start = list(k = 0.5))

  expect_lt(summary(pippo)$parameters[4], 0.05)
  expect_equal(check_fit(pippo), TRUE, ignore_attr = TRUE)
})

test_that("check fit returns true with a significant linear regression", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.9, 0.5, 0.1, 0.8, 0.7)
  pippo <- nls(conc ~ exp(-k * t), start = list(k = 0.5))

  expect_gt(summary(pippo)$parameters[4], 0.05)
  expect_equal(check_fit(pippo), FALSE, ignore_attr = TRUE)
})


# logformerr
test_that("logformerr returns the right output", {
  dframe <- data.frame(alpha = c(1, 2, 3),
                       beta = c(20, 10, 5),
                       err_beta = c(0.5, 0.4, 0.3))
  fdframe <- data.frame(
    log_x = c(log((20 + 10) / 2), log((10 + 5) / 2)),
    log_dx_dt = c(log((20 - 10) / (2 - 1)), log((10 - 5) / (3 - 2))),
    err_x = c(
      (sqrt(0.5^2 + 0.4^2) / 2) / ((10 + 20) / 2),
      (sqrt(0.4^2 + 0.3^2) / 2) / ((10 + 5) / 2)
    ),
    err_log = c(
      (sqrt(0.4^2 + 0.5^2) / (2 - 1)) / ((20 - 10) / (2 - 1)),
      (sqrt(0.4^2 + 0.3^2) / (3 - 2)) / ((10 - 5) / (3 - 2))
    )
  )

  expect_equal(logformerr(dframe), fdframe)
})

# logform
test_that("logform functions as expected", {
  dframe <- data.frame(alpha = c(1, 2, 3), beta = c(20, 10, 5))
  fdframe <- data.frame(
    log_x = c(log((20 + 10) / 2), log((10 + 5) / 2)),
    log_dx_dt = c(log((20 - 10) / (2 - 1)), log((10 - 5) / (3 - 2)))
  )
  expect_equal(logform(dframe), fdframe)
})

# par_est
test_that("par_est estimates correctly the parameter with unweighted data with
          n=1", {
  t <- c(0, 1, 2, 3, 4, 5)
  y1 <- c(128, 64, 32, 16, 8, 4)
  dframe <- data.frame(t, y1)
  expect_no_error(est <- par_est(dframe, 1))
  expect_equal(par_est(dframe, 1), 0.69)

  y_check <- exp(-est * t) * 128
  expect_equal(y1, y_check, tolerance = 0.1)
})

test_that("par_est estimates correctly the parameter with weighted data with
          n=1", {
  t <- c(0, 1, 2, 3, 4, 5)
  y1 <- c(128, 64, 32, 16, 8, 4)
  err <- c(0.1, 0.2, 0.4, 0.1, 0.4, 0.2)
  dframe <- data.frame(t, y1, err)
  expect_no_error(est <- par_est(dframe, 1))
  expect_equal(par_est(dframe, 1), 0.69)

  y_check <- exp(-est * t) * 128
  expect_equal(y1, y_check, tolerance = 0.01)
})

test_that("par_est estimates correctly the parameter with weighted data with
          n=2", {
  t <- c(0, 1, 2, 3, 4, 5)
  y1 <- c(128, 64, 32, 16, 8, 4)
  err <- c(0.1, 0.2, 0.4, 0.1, 0.4, 0.2)
  dframe <- data.frame(t, y1, err)
  expect_no_error(est <- par_est(dframe, 2))
  expect_equal(par_est(dframe, 2), 0.016)

  y_check <- 128 / (1 + 128 * est * t)
  expect_equal(y1, y_check, tolerance = 10)
})


# det_order
test_that("det_order returns a 0-order model with weighted data", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 0.80, 0.69, 0.59, 0.51)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  expect_no_error(res <- det_order(dframe))
  expect_equal(res[[6]], c(Order = 0))
  expect_s3_class(res[[4]], "lm")
  expect_s3_class(res, "ord_res")
})

test_that("det_order returns a 0-order model with unweighted data", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 0.80, 0.69, 0.59, 0.51)
  dframe <- data.frame(t, conc)
  expect_no_error(res <- det_order(dframe))
  expect_equal(res[[6]], c(Order = 0))
  expect_s3_class(res[[4]], "lm")
  expect_s3_class(res, "ord_res")
})

test_that("det_order works with 1-order model with weighted data", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  expect_no_error(res <- det_order(dframe))
  expect_equal(res[[6]], c(Order = 1))
  expect_s3_class(res[[4]], "nls")
  expect_s3_class(res, "ord_res")
})

test_that("det_order works with 1-order model with unweighted data", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)

  dframe <- data.frame(t, conc)
  expect_no_error(res <- det_order(dframe))
  expect_equal(res[[6]], c(Order = 1))
  expect_s3_class(res[[4]], "nls")
  expect_s3_class(res, "ord_res")
})

test_that("det_order works with 2-order model with unweighted data", {
  t <- c(0, 4, 8, 12, 16, 20)

  conc <- 10 / (1 + 10 * 0.5 * t) * (1 + rnorm(length(t), 0, 0.1))

  dframe <- data.frame(t, conc)
  expect_no_error(res <- det_order(dframe))
  expect_equal(res[[6]], c(Order = 2), tolerance = 0.8)
  expect_equal(k_value(res), c(k = 0.5), tolerance = 0.5)
  expect_s3_class(res[[4]], "nls")
  expect_s3_class(res, "ord_res")
})

test_that("det_order works with 2-order model with weighted data", {
  t <- c(0, 4, 8, 12, 16, 20)
  err <- rnorm(length(t), 0, 0.1)
  conc <- 10 / (1 + 10 * 0.5 * t) * (1 + err)

  dframe <- data.frame(t, conc, err)
  expect_no_error(res <- det_order(dframe))
  expect_equal(res[[6]], c(Order = 2), tolerance = 0.8)
  expect_equal(k_value(res), c(k = 0.5), tolerance = 0.5)
  expect_s3_class(res[[4]], "nls")
  expect_s3_class(res, "ord_res")
})

test_that("if no order is found particular output is given", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.51, 0.7, 0.12, 0.07)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  expect_no_error(res <- det_order(dframe))
  expect_equal(res[[6]], c(Order = NaN))
  expect_s3_class(res, "ord_res")
  expect_equal(res[[4]],
    "Ordinary n-order kinetic models do not describe accurately
        the process.",
    ignore_attr = TRUE
  )
})

test_that("det_order returns error with non-dataframe type input", {
  a <- c(1, 2, 3, 4, 5)
  expect_error(det_order(a))
})

test_that("det_order returns error with dataframe with columns numbers other
          than 2 and 3", {
  a <- c(1, 2, 3, 4, 5)
  d1 <- data.frame(a)
  d2 <- data.frame(a, a, a, a)
  expect_error(det_order(d1))
  expect_error(det_order(d2))
})

test_that("det order does not work with a dataframe with less than 3 points", {
  t <- c(0, 4, 8)
  conc <- c(1, 0.51, 0.27)
  err <- c(0.02, 0.05, 0.04)
  dframe <- data.frame(t, conc, err)
  dframe1 <- data.frame(t, conc)
  dframe2 <- data.frame(x = c(1, 2), y = c(3, 4))
  expect_error(res <- det_order(dframe))
  expect_error(res <- det_order(dframe1))
  expect_error(res <- det_order(dframe2))
})

test_that("det_order does not work after duplicates return a datframe with
          less than 3 points", {
  t <- c(0, 4, 8, 12)
  conc <- c(1, 0.51, 0.51, 51)
  err <- c(0.02, 0.05, 0.04, 0.02)
  dframe <- data.frame(t, conc, err)
  expect_error(res <- det_order(dframe))
})


test_that("det_order works with tibble", {
  skip_if_not_installed("tibble")
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- tibble::as_tibble(data.frame(t, conc, err))
  expect_no_error(det_order(dframe))
})
