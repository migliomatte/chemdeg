# Chiquad_red
test_that("Chiquad_red returns error if a non-nls/ord-res object is passed", {
  lr1 <- lm(c(9, 7, 4, 3) ~ c(1, 2, 4, 5))
  pippo <- c(1, 2, 3, 4, 5)
  pippo2 <- data.frame(x = c(1, 2, 3, 4, 5), y = c(1, 2, 4, 5, 25))
  expect_error(chiquad_red(lr1))
  expect_error(chiquad_red(pippo))
  expect_error(chiquad_red(pippo2))
})

test_that("Chiquad_red works with nls weighted object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  regr <- nls(conc ~ exp(-k * t), weights = 1 / err^2, start = list(k = 0.5))
  expect_no_error(chiquad_red(regr))
  expect_equal(round(chiquad_red(regr), 4), round(0.1344415, 4))
})

test_that("Chiquad_red works with ord-res weighted object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_no_error(chiquad_red(res))
  expect_equal(round(chiquad_red(res), 4), round(0.1344415, 4))
})

test_that("Chiquad_red doesn't works with nls unweighted object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  regr <- nls(conc ~ exp(-k * t), start = list(k = 0.5))
  expect_error(chiquad_red(regr))
})

test_that("Chiquad_red doesn't works with ord-res unweighted object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  dframe <- data.frame(t, conc)
  res <- det_order(dframe)
  expect_error(chiquad_red(res))
})

test_that("Chiquad_red returns error with ord-res 0-order", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 0.80, 0.69, 0.59, 0.51)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_error(chiquad_red(res))
})

# AICC
test_that("AICC works with weighted nls object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  regr <- nls(conc ~ exp(-k * t), weights = 1 / err^2, start = list(k = 0.5))
  expect_no_error(AICC(regr))
  expect_equal(signif(AICC(regr), 2), -35, tolerance = 0.1)
})

test_that("AICC works with unweighted nls object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  regr <- nls(conc ~ exp(-k * t), start = list(k = 0.5))
  expect_no_error(AICC(regr))
  expect_equal(signif(AICC(regr), 4), -35.6, tolerance = 0.1)
})

test_that("AICC returns error if a non-nls object is passed", {
  lr1 <- lm(c(9, 7, 4, 3) ~ c(1, 2, 4, 5))
  pippo <- c(1, 2, 3, 4, 5)
  pippo2 <- data.frame(x = c(1, 2, 3, 4, 5), y = c(1, 2, 4, 5, 25))
  expect_error(AICC(lr1))
  expect_error(AICC(pippo))
  expect_error(AICC(pippo2))
})

# goodness_of_fit
test_that("goodness_of_fit returns error if a non-nls/ord-res/lm object is
          passed", {
  pippo <- c(1, 2, 3, 4, 5)
  pippo2 <- data.frame(x = c(1, 2, 3, 4, 5), y = c(1, 2, 4, 5, 25))

  expect_error(goodness_of_fit(pippo))
  expect_error(goodness_of_fit(pippo2))
})

test_that("goodness_of_fit works with a lm object", {
  lr1 <- lm(c(9, 7, 4, 3) ~ c(1, 2, 4, 5))
  expect_no_error(goodness_of_fit(lr1))
})

test_that("goodness_of_fit works with nls weighted object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  regr <- nls(conc ~ exp(-k * t), weights = 1 / err^2, start = list(k = 0.5))
  expect_no_error(goodness_of_fit(regr))
})

test_that("goodness_of_fit works with ord-res weighted object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_no_error(goodness_of_fit(res))
})

test_that("goodness_of_fit works with nls unweighted object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  regr <- nls(conc ~ exp(-k * t), start = list(k = 0.5))
  expect_no_error(goodness_of_fit(regr))
})

test_that("goodness_of_fit works with ord-res unweighted object", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
  dframe <- data.frame(t, conc)
  res <- det_order(dframe)
  expect_no_error(goodness_of_fit(res))
})

test_that("goodness_of_fit returns the summary of fit with ord-res 0-order", {
  t <- c(0, 4, 8, 12, 16, 20)
  conc <- c(1, 0.91, 0.80, 0.69, 0.59, 0.51)
  err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
  dframe <- data.frame(t, conc, err)
  res <- det_order(dframe)
  expect_no_error(goodness_of_fit(res))
})
