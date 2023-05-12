lmp <- function(lmobj) {
  if (methods::is(lmobj, "lm")) {
    f_values <- summary(lmobj)$fstatistic
    pr <- stats::pf(f_values[1], f_values[2], f_values[3], lower.tail = FALSE)
    attributes(pr) <- NULL
    return(pr)
  } else {
    stop("The argument provided is not a 'lm' class object")
  }
}

#' Formula of an n-order model.
#'
#' Given the reaction order \eqn{n} , the function returns the equation
#' corresponding to that particular n^th^-order kinetic model.
#' For \eqn{n\neq 1}:
#' \deqn{y(t)=((n-1)\,k\,t+y_0^{1-n}))^{\frac{1}{n-1}}}
#' for \eqn{n=1}:
#' \deqn{y(t)=y_0\,e^{-k\,t}}
#' @param n reaction order
#'
#' @return A formula object containing the equation of the selected *n*^th^
#' order kinetic model.
#' @importFrom stats as.formula
#' @examples
#' nc <- 2
#' f_gen(nc)
#'
#' f_gen(1)
#' @export
#'
f_gen <- function(n) {
  if (n == 1) {
    form <- y~y0 * exp(-k * t)
  } else {
    n1 <- as.character(n - 1)
    n2 <- as.character(1 / (1 - n))
    f_string <- paste("y~(", n1, "*(k*t+(y0^(-", n1, "))/", n1, "))^",
                      n2, sep = "")
    form <- stats::as.formula(f_string)
  }
  return(form)
}


dup_rem <- function(dfprova) {
  nc <- 1
  while (nc < length(unlist(dfprova[, 1]))) {
    if (dfprova[nc, 2] == dfprova[nc + 1, 2]) {
      nt <- (dfprova[nc, 1] + dfprova[nc + 1, 1]) / 2
      if (length(dfprova[1, ]) == 2) {
        dfprova <- rbind.data.frame(dfprova, c(nt, dfprova[nc, 2]))
      } else {
        nerr <- dfprova[nc, 3] + dfprova[nc + 1, 3]
        dfprova <- rbind.data.frame(dfprova, c(nt, dfprova[nc, 2], nerr))
      }
      dfprova <- dfprova[-c(nc, nc + 1), ]
      dfprova <- dfprova[do.call(base::order, as.list(dfprova)), ]

      nc <- 1
      next
    }
    nc <- nc + 1
  }
  return(dfprova)
}



#' @importFrom methods is
check_fit <- function(obj) {
  if (methods::is(obj, "nls")) {
    if (summary(obj)$parameters[4] < 0.05) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    if (is(obj, "lm")) {
      if (lmp(obj) < 0.05) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }
}


logformerr <- function(dataf) {
  log_dx_dt <- c()
  log_x <- c()
  err_log <- c()
  err_x <- c()
  for (i in 2:length(unlist(dataf[1]))) {
    names(dataf) <- c("t", "x", "err")
    log_dx_dt <- append(
      log_dx_dt,
      log(abs((dataf$x[i] - dataf$x[i - 1]) /
        (dataf$t[i] - dataf$t[i - 1])))
    )
    err_log <- append(
      err_log,
      abs((sqrt(dataf$err[i]^2 + dataf$err[i - 1]^2) /
        (dataf$x[i] - dataf$x[i - 1])))
    )
    log_x <- append(
      log_x,
      log((dataf$x[i] + dataf$x[i - 1]) / 2)
    )
    err_x <- append(
      err_x,
      abs((sqrt(dataf$err[i]^2 + dataf$err[i - 1]^2) / 2) *
        (1 / ((dataf$x[i] + dataf$x[i - 1]) / 2)))
    )
  }

  return(data.frame(
    log_x = log_x,
    log_dx_dt = log_dx_dt,
    err_x = err_x,
    err_log = err_log
  ))
}


logform <- function(data) {
  names(data) <- c("t", "x")
  log_dx_dt <- c()
  log_x <- c()
  for (i in 2:length(unlist(data[1]))) {
    log_dx_dt <- append(
      log_dx_dt,
      log(abs((data$x[i] - data$x[i - 1]) /
        (data$t[i] - data$t[i - 1])))
    )
    log_x <- append(log_x, log((data$x[i] + data$x[i - 1]) / 2))
  }
  df <- data.frame(log_x = log_x, log_dx_dt = log_dx_dt)
  return(df)
}

#' @importFrom stats weighted.mean
par_est <- function(dframe, n) {
  t <- dframe[[1]][-1]
  y <- dframe[[2]][-1]
  y0 <- dframe[[1, 2]]

  if (n == 1) {
    ka <- -log(y / y0) / t
  } else {
    ka <- (y^(1 - n) - (y0^(1 - n))) / ((n - 1) * t) # t!=0
  }

  if (length(dframe[1, ]) == 3) {
    newerr <- abs(dframe[[3]][-1] / (t * y))
    k <- signif(stats::weighted.mean(ka, 1 / newerr^2), 2)
  } else {
    k <- signif(mean(ka), 2)
  }
  return(k)
}

#' Determining reaction order and kinetic formula
#'
#' The functions seeks to determine the reaction order and kinetic rate
#' constant for chemical models that best fit degradation kinetic data.
#' The input of the function is a data-frame organized as follows:
#' \enumerate{
#' \item first columns, time data;
#' \item second columns, concentration data;
#' \item third column (optional, but higly recommended), experimental error}
#'

#'
#' @param dframe a data-frame with 2 or 3 columns, containing time,
#' concentrations, and (optional) error data.
#'
#' @return A `ord_res` object containing in a list the following information:
#' \enumerate{
#' \item the phase space coordinates of transformed data;
#' \item the linear regression performed in the phase space;
#' \item a boolean variable indicating if the estimate of the degradation rate
#'  constant is statistically significant;
#' \item non-linear regression performed using a n^th^-order kinetic model
#' (if n=0 the regression is linear);
#' \item the data-frame given as the input;
#' \item the estimated reaction order.}
#' @export
#'
#' @importFrom stats lm nls confint
#' @importFrom methods is
#' @seealso [results()] to print the results or [goodness_of_fit()] to visualize
#' the major goodness-of-fit measures; [plot_ord()] to plot the regressions in
#' both the phase and conventional spaces;
#' [kin_regr()] to extract the best kinetic model that explain the data and
#' [phase_space()] to extract the linear regression in the phase space.
#' @examples
#' t <- c(0, 4, 8, 12, 16, 20)
#' conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
#' err <- c(0.02, 0.05, 0.04, 0.04, 0.03, 0.02)
#' dframe <- data.frame(t, conc)
#' res <- det_order(dframe)
#'
#' class(res)
#'
#'
#'
#' dframe2 <- data.frame(t, conc, err)
#' res2 <- det_order(dframe2)
#'
#' res2[[5]] == dframe2
#'
det_order <- function(dframe) {
  if (is.data.frame(dframe)) {
    if (methods::is(dframe, "tbl_df")) {
      dframe <- as.data.frame(dframe)
    }
    tdframe <- dup_rem(dframe)


    if (length(unlist(dframe[1])) <= 3) {
      stop("Error: only 3 or less points in dataset. Minimum requirement is 4")
    }

    y0 <- dframe[1, 2]


    if (length(dframe[1, ]) == 2) {
      nc <- 2
      tdframe <- cbind.data.frame(tdframe[1], tdframe[2] / y0)
      ldframe <- logform(tdframe)

      lmdframe <- stats::lm(log_dx_dt ~ log_x, data = ldframe)
    } else {
      if (length(dframe[1, ]) == 3) {
        nc <- 3
        tdframe <- cbind.data.frame(tdframe[1], tdframe[2:3] / y0)
        ldframe <- logformerr(tdframe)

        lmdframe <- stats::lm(log_dx_dt ~ log_x,
          weights = 1 / ldframe$err_log^2,
          data = ldframe
        )
      } else {
        stop("Error: the number of columns provided does not match the
             requirements. Only dataframe with 2 columns (time,data)
             and 3 columns (time,data,error) are accepted")
      }
    }

    if (lmp(lmdframe) < 0.05) {
      uu <- stats::confint(lmdframe)[2, 1]
      ll <- stats::confint(lmdframe)[2, 2]
      if (ceiling(uu) > floor(ll)) {
        n <- round(lmdframe$coefficients[2], 1)
      } else {
        n <- round(lmdframe$coefficients[2], 0)
      }
      if (n < 0) {
        n <- 0
      }
      fcin <- f_gen(n)


      names(n) <- "Order"
      message("Reaction order estimated: ", n)

      if (nc == 2) {
        nlsdframe <- stats::nls(fcin,
          data = list(
            t = dframe[[1]],
            y = dframe[[2]],
            y0 = y0
          ),
          start = list(
            k = par_est(dframe, n)
          )
        )
      }
      if (nc == 3) {
        err <- 1 / dframe[[3]]^2
        nlsdframe <- stats::nls(fcin,
          data = list(
            t = dframe[[1]],
            y = dframe[[2]],
            y0 = y0,
            err = err
          ),
          weights = err,
          start = list(
            k = par_est(dframe, n)
          )
        )
      }

      st_finalfit <- check_fit(nlsdframe)
      rt <- list(ldframe, lmdframe, st_finalfit, nlsdframe, dframe, n)
      class(rt) <- "ord_res"
      return(rt)
    } else {
      if (summary(lmdframe)$coefficients[1, 4] < 0.05) {
        n <- 0
        names(n) <- "Order"
        message("Reaction order estimated: ", n)

        if (nc == 2) {
          lm0dframe <- stats::lm(y ~ t,
            data = list(
              t = dframe[[1]],
              y = dframe[[2]]
            )
          )
        }

        if (nc == 3) {
          err <- 1 / dframe[[3]]^2
          lm0dframe <- stats::lm(y ~ t,
            data = list(
              t = dframe[[1]],
              y = dframe[[2]],
              err = err
            ),
            weights = err
          )
        }
        st_finalfit <- check_fit(lm0dframe)
        rt <- list(ldframe, lmdframe, st_finalfit, lm0dframe, dframe, n)
        class(rt) <- "ord_res"
        return(rt)
      } else {
        st_finalfit <- FALSE
        lm0dframe <- "Ordinary n-order kinetic models do not describe accurately
        the process."
        n <- NaN
        names(n) <- "Order"

        rt <- list(ldframe, lmdframe, st_finalfit, lm0dframe, dframe, n)
        class(rt) <- "ord_res"
        return(rt)
      }
    }
  } else {
    stop("Error: the argument provided is not a data frame type.")
  }
}
