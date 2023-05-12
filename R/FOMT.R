#' First-Order Multi-Target model regression
#'
#' The function performs a non-linear regression using the first-order
#' multi-target model. The model equation is:
#' \deqn{\frac{S}{S_0}=1-(1-e^{-k\,t})^m}
#' where \eqn{S/S_0} is the fraction of surviving molecules, \eqn{k} is the
#' average number of hits per time unit, \eqn{m} is the number of hits required
#'  to degrade the molecule, and  \eqn{t} is time.
#'
#' @param dtframe A data-frame containing 2 or 3 columns: time, normalized
#' concentration and error (optional), respectively
#'
#' @return Returns the results of the regression as a
#' [nls][stats::nls()] object.
#'
#' @details
#' The FOMT model has been proposed as an alternative to the Weibull equation
#' that is commonly used when the time-dependent behavior of the data
#' significantly deviates from that predicted by standard chemical models.
#'
#' @seealso [FOMTm()][FOMTm()], [par_est_FOMT()]
#' @export
#' @importFrom stats nls
#'
#' @examples
#' t <- c(0, 4, 8, 12, 16, 20)
#' conc <- c(1, 0.98, 0.99, 0.67, 0.12, 0.03)
#' err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
#' dframe <- data.frame(t, conc, err)
#' FOMT <- FOMT(dframe)
#' plot(dframe[[1]], dframe[[2]])
#' arrows(dframe[[1]], dframe[[2]] + dframe[[3]],
#'   dframe[[1]], dframe[[2]] - dframe[[3]],
#'   length = 0
#' )
#' newt <- seq(0, 21, by = 0.1)
#' lines(newt, predict(FOMT, newdata = list(t = newt)))
#'
#' dframe1 <- data.frame(t, conc)
#' FOMT1 <- FOMT(dframe1)
#' plot(dframe1[[1]], dframe1[[2]])
#' lines(newt, predict(FOMT1, newdata = list(t = newt)))
#' summary(FOMT)
#' summary(FOMT1)
FOMT <- function(dtframe) {
  t <- dtframe[[1]]
  y <- dtframe[[2]]

  if (length(dtframe[1, ]) == 2) {
    FOMT <- stats::nls(y ~ 1 - (1 - exp(-k * t))^n,
      data = list(y = y, t = t),
      start = par_est_FOMT(t, y)
    )
  } else {
    if (length(dtframe[1, ]) == 3) {
      err <- dtframe[[3]]
      FOMT <- stats::nls(y ~ 1 - (1 - exp(-k * t))^n,
        data = list(y = y, t = t),
        start = par_est_FOMT(t, y),
        weights = 1 / err^2
      )
    } else {
      stop("Unexpected number of columns")
    }
  }
  return(FOMT)
}


#' First-Order Multi-Target model
#'
#' Call the function to return the formula of the Single-Hit Multi-Target
#' model (FOMT):
#' \deqn{1-(1-e^{-k\,t})^n}
#'
#' @param k average number of hits per time unit
#' @param t time
#' @param m minimum number of hits required to degrade the
#' molecule
#'
#' @return Returns calculated values using the formula of FOMT model
#'
#' It can be used inside the function [nls][stats::nls()] as the RHS of the
#' formula.
#' @export
#'
#' @seealso [FOMT()], [par_est_FOMT()], [stats::nls()]
#'
#' @examples
#' t <- seq(0, 100, by = 1)
#' k <- 0.1
#' n <- 200
#' y <- FOMTm(t, k, n)
#' plot(t, y, type = "l")
FOMTm <- function(t, k, m) {
  1 - (1 - exp(-k * t))^m
}

#' First-Order Multi-Target parameter starting values
#'
#' `par_est_FOMT` estimates the starting values of the parameters of the
#' first-order multi-target model from a data-set.
#'
#' @param x A data-frame with time and concentration in the first and second
#' columns,respectively. Alternatively, it could be an array of time and y an
#' array of concentrations.
#' @param y Optional, an array of concentrations. To be inserted only if x is an
#' array.
#'
#' @return The function returns an array with the suggested initial values of
#' parameters.
#' @export
#' @importFrom stats lm coef
#' @seealso [FOMT()], [FOMTm()]
#' @examples
#' t <- seq(0, 30, by = 6)
#' k <- 0.3
#' n <- 40
#' set.seed(100)
#' y <- FOMTm(t, k, n) * (1 + rnorm(length(t), 0, 0.05))
#'
#' nlsFOMT <- nls(y ~ FOMTm(t, k, n),
#'   data = list(y = y, t = t),
#'   start = par_est_FOMT(t, y)
#' )
#' summary(nlsFOMT)
par_est_FOMT <- function(x, y = NULL) {
  if (is.data.frame(x)) {
    t <- x[[1]]
    y <- x[[2]]
  } else {
    t <- x
  }
  t <- t[y < 0.2]
  y <- y[y < 0.2]

  coefs <- stats::coef(
    stats::lm(I(log(y)) ~ t)
  )
  coefs[1] <- exp(coefs[1])
  coefs[2] <- -coefs[2]
  attributes(coefs) <- NULL

  return(c(k = coefs[2], n = coefs[1]))
}
