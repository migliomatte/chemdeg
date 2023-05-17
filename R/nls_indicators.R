#' Reduced chi-squared
#'
#' Function that returns the reduced chi-squared (\eqn{\chi^2_{red}=\chi^2/df},
#' where \eqn{df} are the degrees of freedom) value for
#' a non-linear regression model (`nls` object). Reduced-chi squared is a
#' goodness-of-fit measure. Values close to 1 indicates a good fit, while
#' values \eqn{>>1} indicate poor fit and values \eqn{<1} indicate
#' over-fitting.
#' The function is calculated only with non-linear regression weighted on
#' experimental error.
#'
#' @param fit `nls` object with weighted fit
#'
#' @return Returns the reduced chi-squared value
#' @export
#' @seealso [stats::dchisq()] for chi-squared distribution; [stats::AIC()],
#' [stats::BIC()], [stats::sigma()] (for RMSE), [AICC()] for other
#' goodness-of-fit
#' indicators. [goodness_of_fit()]
#'
#' @importFrom methods is
#' @importFrom stats residuals
#'
#' @references
#' Philip R. Bevington,
#' D. Keith Robinson,
#' J. Morris Blair,
#' A. John Mallinckrodt,
#' Susan McKay (1993).
#' *Data Reduction and Error Analysis for the Physical Sciences*
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(1.2, 3.9, 8.6, 17.4, 26)
#' er <- c(0.5, 0.8, 0.5, 1.9, 1.2)
#' fit1 <- nls(y ~ k * x^2,
#'   data = list(x = x, y = y),
#'   start = list(k = 1),
#'   weights = 1 / er^2
#' )
#' chiquad_red(fit1)
#'
#' fit2 <- nls(y ~ k * x^3,
#'   data = list(x = x, y = y),
#'   start = list(k = 1),
#'   weights = 1 / er^2
#' )
#' chiquad_red(fit2)
chiquad_red <- function(fit) {
  if (methods::is(fit, "ord_res")) {
    fit <- fit[[4]]
  }
  if (!methods::is(fit, "nls")) {
    stop("Invalid argument: only 'nls' object accepted")
  }
  x <- length(fit$m$lhs())
  if (is.null(fit$weights)) {
    stop("Non-weighted fit was provided. Chi-squared statistic is computed on
         weighted data.")
  }
  se <- fit$weights^-0.5
  residui <- stats::residuals(fit)[1:x]
  VV <- fit$m$gradient()
  Q1 <- qr.Q(qr(VV))
  H <- Q1 %*% t(Q1)
  hh <- diag(H)
  seresidui <- se * sqrt(1 - hh)
  residuistandard <- residui / seresidui
  chiquad <- sum(residuistandard^2)
  df <- x - 2
  chiquadred <- chiquad / df

  return(chi_red = chiquadred)
}

#' Akaike Information Criterion With Correction
#'
#' The function calculates the Akaike Information Criterion with correction
#' for small samples size.
#'
#' @param fit a 'nls'-object
#' @export
#' @return Returns the AICc value
#'
#' @details
#' When the sample size is small, there is a substantial probability that AIC
#' (see [stats::AIC()] for more details)
#' will select models that have too many parameters, i.e. that AIC will
#' overfit. AICc is AIC with a correction for small sample sizes.
#'
#' The AICc is computed as follows:
#' \deqn{AICc=AIC+\frac{2\,k\,(k+1)}{n-k-1}}
#' where n denotes the sample size and k denotes the number of parameters. Thus
#' , AICc is essentially AIC with an extra penalty term for the number of
#' parameters. Note that as \eqn{n\rightarrow \infty}, the extra penalty term
#' converges to 0, and thus AICc converges to AIC.
#'
#' @seealso [stats::AIC()] for uncorrected AIC, [stats::BIC()],
#' [stats::sigma()] ,[chiquad_red()] for other goodness of fit indicators.
#' [goodness_of_fit()]
#' @importFrom stats AIC
#' @examples
#' t <- seq(0, 10, 1)
#' y <- 1 / (0.5 * exp(t) + 1) + stats::rnorm(length(t), 0, 0.05)
#'
#' fit <- nls(y ~ 1 / (k * exp(t) + 1),
#'   data = list(t = t, y = y),
#'   start = list(k = 0.2)
#' )
#' AICC(fit)
#'
AICC <- function(fit) {
  if (methods::is(fit, "nls")) {
    k <- length(fit$m$getAllPars())
    n <- length(fit$m$lhs())
    aicc <- stats::AIC(fit) + 2 * k * (k + 1) / (n - k - 1)
    return(aicc)
  } else {
    stop("Invalid argument: only 'nls' object accepted")
  }
}

#' Goodness-of-fit, non-linear regression
#'
#' Function that returns the following goodness-of-fit statistics for
#' non-linear regression: AIC, AICc, BIC, RMSE and reduced Chi-squared.
#'
#' The function returns the values of AIC, AICC, BIC, RMSE and reduced
#' chi-squared (\eqn{\chi^2_{red}}) for `nls` objects. If a linear model
#' object is passed, the function returns its [summary][stats::summary.lm()].
#'
#' Given an `ord_res` object (output of the function [det_order()]),
#' the function returns one of the results
#' above depending on the model chosen to explain the data.
#'
#' Because the [chiquad_red()] function returns the value only with weighted
#' data, the \eqn{\chi^2_{red}} will be returned only with weighted
#' regressions.
#'
#' @param fit a `nls`, `lm` or ord_res object
#'
#' @return It returns a table with the values of AIC, AICc, BIC, RSME and
#' reduced Chi squared. Single goodness-of-fit measures can  be obtained as
#' follows:
#' \enumerate{
#' \item call standard R functions [stats::AIC()], [stats::BIC()],
#' [stats::sigma()] for AIC, BIC and RMSE, respectively;
#' \item call `chemdeg` functions [AICC()] and [chiquad_red()] for AICc and
#' reduced chi-squared, respectively.
#' }
#' @export
#' @seealso [stats::AIC()], [AICC()], [stats::BIC()], [stats::sigma()],
#' [chiquad_red()]
#' @importFrom stats AIC BIC sigma
#' @importFrom methods is
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(1.2, 3.9, 8.6, 17.4, 26)
#' er <- c(0.5, 0.8, 0.5, 1.9, 1.2)
#' fit1 <- nls(y ~ k * x^2,
#'   data = list(x = x, y = y), start = list(k = 1),
#'   weights = 1 / er^2
#' )
#' goodness_of_fit(fit1)
goodness_of_fit <- function(fit) {
  if (methods::is(fit, "ord_res")) {
    fit <- fit[[4]]
  }
  if (methods::is(fit, "nls")) {
    aic <- stats::AIC(fit)
    aicc <- AICC(fit)
    bic <- stats::BIC(fit)
    oldoptions<-options()
    on.exit(options(oldoptions))
    options(show.error.messages = FALSE)
    chi_red <- try(chiquad_red(fit))
    if (!methods::is(chi_red, "numeric")) {
      chi_red <- NA
    }
    rmse <- stats::sigma(fit)

    tb <- matrix(c(aic, aicc, bic, rmse, chi_red), ncol = 1)
    colnames(tb) <- "Value"
    row.names(tb) <- c("AIC:", "AICc:", "BIC:", "RMSE:", "Chi-sq_red:")
    return(tb)
  } else {
    if (methods::is(fit, "lm")) {
      return(summary(fit))
    } else {
      stop("Invalid object class. Only \"lm\", \"nls\" and \"ord-res\"
           class objects are accepted.")
    }
  }
}
