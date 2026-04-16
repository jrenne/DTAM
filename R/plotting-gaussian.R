#' Plot Gaussian iso-density contours
#'
#' Draws contour lines associated with a bivariate Gaussian distribution.
#'
#' @param Ew Mean vector of length two.
#' @param Vw `2 x 2` covariance matrix.
#' @param n Grid size used for contour evaluation.
#' @param min_w,max_w Optional plotting bounds. If omitted, the routine uses a
#'   multiple of the marginal standard deviations.
#' @param prob_levels Probability-content levels used to determine contour
#'   radii.
#' @param n_std Number of standard deviations used when bounds are not supplied.
#' @param maintitle Plot title.
#'
#' @return Invisibly returns a list containing the plotting grid and contour
#'   levels.
#'
#' @examples
#' plot_isodensity_gaussian(Ew = c(0, 0),
#'                          Vw = matrix(c(1, 0.5, 0.5, 2), 2, 2),
#'                          maintitle = "")
#'
#' @export
plot_isodensity_gaussian <- function(Ew, Vw, n = 200,
                                     min_w = NA_real_, max_w = NA_real_,
                                     prob_levels = c(0.50, 0.80, 0.95, 0.99),
                                     n_std = 4,
                                     maintitle = "") {
  stopifnot(length(Ew) == 2, all(dim(Vw) == c(2, 2)))

  sds <- sqrt(diag(Vw))

  # Default axis range (if not supplied)
  if (!is.finite(min_w)) min_w <- Ew[1] - n_std * sds[1]
  if (!is.finite(max_w)) max_w <- Ew[1] + n_std * sds[1]

  # Grid for density
  x1 <- seq(Ew[1] - n_std * sds[1], Ew[1] + n_std * sds[1], length.out = n)
  x2 <- seq(Ew[2] - n_std * sds[2], Ew[2] + n_std * sds[2], length.out = n)

  X1 <- matrix(rep(x1, each = n), nrow = n)
  X2 <- matrix(rep(x2, times = n), nrow = n)

  detV <- det(Vw)
  invV <- solve(Vw)

  dX1 <- X1 - Ew[1]
  dX2 <- X2 - Ew[2]

  # Quadratic form Q = (w-Ew)' V^{-1} (w-Ew)
  Q <- invV[1,1]*dX1^2 + 2*invV[1,2]*dX1*dX2 + invV[2,2]*dX2^2

  # Density
  logdens <- -0.5 * (2 * log(2*pi) + log(detV) + Q)
  dens <- exp(logdens)

  # Convert probability-mass levels to density levels
  # In 2D Gaussian: Q ~ chi^2(df=2) on ellipses
  Q_levels <- qchisq(prob_levels, df = 2)
  dens_levels <- exp(-0.5 * Q_levels) / (2*pi*sqrt(detV))

  contour(x1, x2, dens,
          las = 1,
          xlim = c(min_w, max_w),
          ylim = c(min_w, max_w),
          levels = dens_levels,
          labels = paste0(100 * prob_levels, "%"),
          drawlabels = TRUE,
          xlab = expression(r[t]),
          ylab = expression(x[t]),
          main = maintitle,
          labcex = 0.8,
          vfont = c("sans serif", "bold"),
          lwd = 2, col = "grey")

  grid()

  # Mark unconditional mean
  points(Ew[1], Ew[2], pch = 3, lwd = 2)
  text(Ew[1], Ew[2], labels = expression(E[t](w[T])), pos = 4)

  invisible(list(x1 = x1, x2 = x2, dens = dens,
                 dens_levels = dens_levels, prob_levels = prob_levels))
}
