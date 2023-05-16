#' Phase space, linearized model
#'
#' Given an `ord_res` object, this function returns the linearized model
#' that best fits the data in the phase space.
#' `ord_res` object can be obtained using the function [det_order()].
#'
#' @param x an `ord_res` object
#'
#' @return Returns a `lm` class object.
#' @seealso [det_order()], [kin_regr()], [results()], [stats::lm()]
#' @export
#' @examples
#' t <- c(0, 4, 8, 12, 16, 20)
#' conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
#' dframe <- data.frame(t, conc)
#' res <- det_order(dframe)
#'
#' phase_space(res)
phase_space <- function(x) {
  stopifnot(methods::is(x,"ord_res"))
  return(x[[2]])
}

#' Degradation kinetics
#'
#' Returns from an `ord_res` object either the linear or the non-linear
#' regression of the degradation kinetics data.
#'
#' After the
#' analysis in the phase space for the determination of the reaction order,
#' [det_order()] performs either a linear or a non-linear
#' regression of the kinetic data, depending on whether the reaction order is
#' n=0 or n>0,
#' respectively. To access the regression object call `kin_degr`.
#'
#' @param x an `ord_res` object
#'
#' @return Returns either an `nls` or `lm` object based on the regression
#' performed by the function [det_order()].
#' @export
#' @seealso [det_order()], [phase_space()], [results()], [stats::lm()]
#' @examples
#' t <- c(0, 4, 8, 12, 16, 20)
#' conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
#' dframe <- data.frame(t, conc)
#' res <- det_order(dframe)
#'
#' kin_regr(res)
kin_regr <- function(x) {
  stopifnot(methods::is(x,"ord_res"))
  return(x[[4]])
}

#' Summary of 'ord_res' object
#'
#' Returns the results of the analyses performed by [det_order()] function.
#'
#' The function prints:
#' \enumerate{
#'    \item{the linear regression performed in the phase space, together with
#'    the estimated *n* value and its 95% confidence interval}
#'    \item{a brief conclusion on the results obtained in the phase space
#'    stating which reaction order should be preferred}
#'    \item{the (non-)linear regression performed with parameters associated
#'    statistics. If a non-linear regression has been performed, the most
#'    common goodness-of-fit measures calculated with [goodness_of_fit()]
#'    are printed}
#'    }
#'
#' @param object an 'ord_res' object
#'
#' @return It prints a summary of the analysis in the phase space, the reaction
#' order, and the regression results.
#' @export
#' @importFrom stats coefficients confint
#' @import MASS
#' @seealso [det_order()], [kin_regr()], [phase_space()]
#' @examples
#' t <- c(0, 4, 8, 12, 16, 20)
#' conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
#' err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
#' dframe <- data.frame(t, conc, err)
#' res <- det_order(dframe)
#'
#' results(res)
results <- function(object) {
  cat("\nLinear regression in the phase space:", "\n")
  cat(
    "log(dx/dt)=", round(object[[2]]$coefficients[[2]], 2), "log(x) + (",
    round(object[[2]]$coefficients[[1]], 2), ")\n\n"
  )
  cat("Estimate of n:\n\n")

  ic <- summary(object[[2]])$coefficients[2, ]

  print(ic,
    digits = getOption("digits"), quote = FALSE,
    na.print = "", zero.print = "0",
    right = is.numeric(ic) || is.complex(ic),
    justify = "none"
  )

  cat("\nConfidence interval of n:", "\n")

  ic <- stats::confint(object[[2]])[2, ]

  print(ic,
    digits = getOption("digits"), quote = FALSE,
    na.print = "", zero.print = "0",
    right = is.numeric(ic) || is.complex(ic),
    justify = "none"
  )

  cat("\n")
  if (is.nan(object[[6]])) {
    cat(
      "The process is not appropriately described by a standard n-order",
      "kinetic model. Alternative models should be evaluated.\n"
    )
    proceed <- FALSE
  } else {
    if (round(object[[6]]) == object[[6]] && object[[6]] != 0) {
      cat("Statistical analysis indicates that an order ", object[[6]],
        " degradation kineitc model is likely to describe the data.\n",
        "The null hypothesis H0:\n",
        "\"The process is described by an order ", object[[6]],
        "kinetic model\"\n cannot be rejected.",
        sep = ""
      )
      proceed <- TRUE

      if (!object[[3]]) {
        cat(
          "\n Estimates of model parameters are not significant.",
          "Alternative models should be evaluated"
        )
        return()
      }
    } else {
      if (object[[6]] == 0) {
        proceed <- FALSE
        cat("Estimate of the intercept is significant; but the estimate of",
          " slope is not.\nThe data are likely to be described by an 0-order",
          " kinetic model.",
          sep = ""
        )

        if (!object[[3]]) {
          cat(
            "\n Estimates of model parameters are not significant.",
            "Alternative models should be evaluated."
          )
          return()
        }

        cat(
          "\n\nLinear  regression was performed (0-order kinetics).\n ",
          "\nEstimated k value: \n"
        )
        print(summary(object[[4]]))
        cat("\nConfidence interval of k: \n")
        ic <- -stats::confint(object[[4]])[2, ]
        print(ic,
          digits = getOption("digits"), quote = FALSE,
          na.print = "", zero.print = "0",
          right = is.numeric(ic) || is.complex(ic),
          justify = "none"
        )
      } else {
        cat(
          "The confidence interval of n excludes integer-valued orders.\n",
          "The best estimate n=", object[[6]], "was used to perform the",
          "regression.\n \n"
        )
        proceed <- TRUE
      }
    }
  }

  if (proceed) {
    cat(
      "\n\nNon-linear least squares regression was performed with an order ",
      object[[6]], " kinetic model:\n \n Estimate of k: \n"
    )
    ic <- summary(object[[4]])$parameters
    print(ic,
      digits = getOption("digits"), quote = FALSE,
      na.print = "", zero.print = "0",
      right = is.numeric(ic) || is.complex(ic),
      justify = "none"
    )

    ci2 <- stats::confint(object[[4]])
    cat("Confidence interval of k: \n")
    print(ci2,
      digits = getOption("digits"), quote = FALSE,
      na.print = "", zero.print = "0",
      right = is.numeric(ci2) || is.complex(ci2),
      justify = "none"
    )


    cat("\nGoodness-of-fit:\n")
    tb <- goodness_of_fit(object[[4]])
    print(tb)
    if (is.na(tb[5, 1])) {
      cat(
        "NB: Reduced Chi-squared is not calculated with",
        "unweighted data\n"
      )
    }
  }

  cat("-----------------------------------------------------\n\n")
}


k_value <- function(ordres) {
  if (is(ordres[[4]], "nls")) {
    return(stats::coefficients(ordres[[4]]))
  } else {
    if (is(ordres[[4]], "lm")) {
      return(abs(stats::coefficients(ordres[[4]])))
    } else {
      message(
        "Kinetic model was not found in object. \nProbably the phase",
        "space regression was not significant."
      )
    }
  }
}


n_value <- function(ordres) {
  stopifnot(methods::is(ordres,"ord_res"))
  return(ordres[[6]])
}
