#' Plot of phase space and kinetic curve
#'
#' The function plots the results obtained from a [det_order()] function.
#' Two plots are shown: one representing the transformed data in the phase
#' space and the other the kinetic data in the conventional space along with
#' their regression curves.
#'
#' @param ord_res an 'ord_res' object
#'
#' @return Two plots. The first representing the transformed data in the phase
#' space and the other the kinetic data in the conventional space along with
#' their regression curves.
#' Black line represent the best regression curve, whereas green lines show the
#' fits with the reaction order chosen.
#' @export
#'
#' @importFrom graphics par arrows plot abline lines
#' @importFrom stats predict
#' @importFrom methods is
#'
#' @examples
#' t <- c(0, 4, 8, 12, 16, 20)
#' conc <- c(1, 0.51, 0.24, 0.12, 0.07, 0.02)
#' dframe <- data.frame(t, conc)
#' res <- det_order(dframe)
#'
#' plot_ord(res)
plot_ord <- function(ord_res) {
  oldpar<-par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (methods::is(ord_res, "ord_res")) {
    dflplot <- as.data.frame(ord_res[[1]])

    dfline <- data.frame(
      c(ord_res[[2]]$coefficients[1]),
      c(ord_res[[2]]$coefficients[2])
    )
    names(dfline) <- c("k", "n")

    if (is.nan(ord_res[[6]])) {
      graphics::plot(dflplot$log_x, dflplot$log_dx_dt,
        main = "Phase space",
        xlab = "log(x)", ylab = "log(|dx/dy|)"
      )
      graphics::abline(a = dfline$n, b = dfline$k)
      message <- "The process is not well described by a n-order kinetic model"
      return(list(message))
    }
    dfplot <- as.data.frame(ord_res[[5]])
    if (length(ord_res[[1]][1, ]) == 2) {
      names(dfplot) <- c("x", "y")
    } else {
      names(dfplot) <- c("x", "y", "err")
    }
    xn <- seq(0, max(dfplot$x) * 1.1, by = max(dfplot$x) / 100)
    pn <- stats::predict(ord_res[[4]], newdata = list(t = xn, y0 = dfplot$y[1]))
    dfnls <- data.frame(xn, pn)

    if (ord_res[[6]] == 0) {
      dfnlord <- data.frame(
        c(ord_res[[6]]),
        c(summary(ord_res[[4]])$coefficients[2])
      )
    } else {
      dfnlord <- data.frame(
        c(ord_res[[6]]),
        c(summary(ord_res[[4]])$parameters[1])
      )
    }
    names(dfnlord) <- c("n", "k")
    mm <- min(dflplot$log_x)
    MM <- min(dflplot$log_dx_dt)

    mm1 <- max(dfplot$x)
    MM1 <- max(dfplot$y)

    graphics::par(pty = "s", ask = TRUE)
    graphics::plot(dflplot$log_x, dflplot$log_dx_dt,
      xlab = "log(x)", ylab = "log(|dx/dy|)", main = "Phase space",
      ylim = c(min(dflplot$log_dx_dt), 0), xlim = c(min(dflplot$log_x), 0)
    )
    graphics::abline(a = dfline$k, b = dfline$n)
    graphics::abline(a = log(abs(dfnlord$k)), b = dfnlord$n, col = "red")

    if (length(ord_res[[1]][1, ]) == 4) {
      graphics::arrows(dflplot$log_x, dflplot$log_dx_dt + dflplot$err_log,
        dflplot$log_x, dflplot$log_dx_dt - dflplot$err_log,
        length = 0
      )
      graphics::arrows(dflplot$log_x + dflplot$err_x, dflplot$log_dx_dt,
        dflplot$log_x - dflplot$err_x, dflplot$log_dx_dt,
        length = 0
      )
    }

    graphics::plot(dfplot$x, dfplot$y,
      xlab = "time", ylab = "Concentration", main = "Kinetic degradation"
    )
    graphics::lines(dfnls$xn, dfnls$pn, col = "red")

    if (length(ord_res[[1]][1, ]) == 4) {
      graphics::arrows(dfplot$x, dfplot$y + dfplot$err,
        dfplot$x, dfplot$y - dfplot$err,
        length = 0
      )
    }

    if (!ord_res[[3]]) {
      message("NB: parameters of the regression are not significant!")
    }
  } else {
    stop("Error: the argument is not an object from det_order function")
  }
}
