make_sigma_points <- function(mean, covariance, alpha, beta, kappa) {
  mean <- as.numeric(mean)
  covariance <- as.matrix(covariance)
  n <- length(mean)
  lambda <- alpha^2 * (n + kappa) - n
  scale <- n + lambda
  if (!is.finite(scale) || scale <= 0) {
    stop("alpha and kappa must imply n + lambda > 0.")
  }

  root <- t(chol(scale * covariance))
  points <- matrix(mean, n, 2L * n + 1L)
  points[, 1L + seq_len(n)] <- points[, 1L + seq_len(n), drop = FALSE] + root
  points[, 1L + n + seq_len(n)] <-
    points[, 1L + n + seq_len(n), drop = FALSE] - root

  weights.mean <- c(lambda / scale, rep(1 / (2 * scale), 2L * n))
  weights.cov <- weights.mean
  weights.cov[1] <- weights.cov[1] + 1 - alpha^2 + beta

  list(points = points,
       weights.mean = weights.mean,
       weights.cov = weights.cov)
}

#' Unscented Kalman filter
#'
#' Applies the scaled unscented Kalman filter to a nonlinear state-space model
#' with additive Gaussian transition and measurement noise:
#' \deqn{x_t = f(x_{t-1},t) + \varepsilon_t,\quad
#' \varepsilon_t \sim N(0,Q),}
#' \deqn{y_t = h(x_t,t) + \eta_t,\quad \eta_t \sim N(0,R).}
#'
#' @param Y Matrix of observations, with one row per date.
#' @param transition Function `transition(x, t)` returning the conditional mean
#'   of the next state.
#' @param measurement Function `measurement(x, t)` returning the conditional
#'   mean of the observation.
#' @param Q State-innovation covariance matrix.
#' @param R Measurement-innovation covariance matrix.
#' @param x0 State mean immediately before the first transition.
#' @param P0 Covariance matrix immediately before the first transition.
#' @param alpha,beta,kappa Parameters of the scaled unscented transform. For a
#'   Gaussian prior, `beta = 2` is conventional.
#'
#' @return A list containing filtered and predicted state means
#'   (`filtered.mean`, `predicted.mean`), the corresponding covariance arrays
#'   (`filtered.cov`, `predicted.cov`), innovations and their covariance arrays
#'   (`innovation`, `innovation.cov`), and the total and date-specific Gaussian
#'   quasi-log-likelihood (`loglik`, `loglik.vector`).
#'
#' @examples
#' # Nonlinear measurement equation: y_t = x_t + 0.1 x_t^2 + eta_t.
#' transition <- function(x, t) 0.8 * x
#' measurement <- function(x, t) x + 0.1 * x^2
#'
#' set.seed(123)
#' TT <- 30
#' x <- numeric(TT)
#' y <- numeric(TT)
#' previous <- 0
#' for (t in seq_len(TT)) {
#'   x[t] <- transition(previous, t) + rnorm(1, sd = 0.2)
#'   y[t] <- measurement(x[t], t) + rnorm(1, sd = 0.1)
#'   previous <- x[t]
#' }
#'
#' fit <- UKF(matrix(y), transition, measurement,
#'            Q = matrix(0.2^2), R = matrix(0.1^2),
#'            x0 = 0, P0 = matrix(0.2))
#' head(fit$filtered.mean)
#'
#' @export
UKF <- function(Y, transition, measurement, Q, R, x0, P0,
                alpha = 1e-3, beta = 2, kappa = 0) {
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  R <- as.matrix(R)
  x0 <- as.numeric(x0)
  P0 <- as.matrix(P0)

  TT <- nrow(Y)
  m <- ncol(Y)
  n <- length(x0)
  if (TT < 1L) {
    stop("Y must contain at least one row.")
  }
  if (anyNA(Y)) {
    stop("UKF does not currently support missing observations.")
  }
  if (any(dim(Q) != c(n, n)) || any(dim(P0) != c(n, n))) {
    stop("Q and P0 must be square matrices matching the state dimension.")
  }
  if (any(dim(R) != c(m, m))) {
    stop("R must be a square matrix matching the observation dimension.")
  }

  filtered.mean <- matrix(NA_real_, TT, n)
  predicted.mean <- matrix(NA_real_, TT, n)
  innovation <- matrix(NA_real_, TT, m)
  filtered.cov <- array(NA_real_, c(n, n, TT))
  predicted.cov <- array(NA_real_, c(n, n, TT))
  innovation.cov <- array(NA_real_, c(m, m, TT))
  loglik.vector <- numeric(TT)

  x.previous <- x0
  P.previous <- P0

  for (t in seq_len(TT)) {
    sigma <- make_sigma_points(x.previous, P.previous, alpha, beta, kappa)
    state.points <- matrix(vapply(
      seq_len(ncol(sigma$points)),
      function(i) as.numeric(transition(sigma$points[, i], t)),
      numeric(n)
    ), nrow = n)
    x.pred <- c(state.points %*% sigma$weights.mean)
    centered.state <- state.points - x.pred
    P.pred <- Q + (centered.state * rep(sigma$weights.cov, each = n)) %*%
      t(centered.state)
    P.pred <- (P.pred + t(P.pred)) / 2

    predicted.sigma <- make_sigma_points(x.pred, P.pred, alpha, beta, kappa)
    measurement.points <- matrix(vapply(
      seq_len(ncol(predicted.sigma$points)),
      function(i) as.numeric(measurement(predicted.sigma$points[, i], t)),
      numeric(m)
    ), nrow = m)
    y.pred <- c(measurement.points %*% predicted.sigma$weights.mean)
    centered.measurement <- measurement.points - y.pred
    centered.predicted.state <- predicted.sigma$points - x.pred
    S <- R +
      (centered.measurement * rep(predicted.sigma$weights.cov, each = m)) %*%
      t(centered.measurement)
    S <- (S + t(S)) / 2
    cross.cov <-
      (centered.predicted.state * rep(predicted.sigma$weights.cov, each = n)) %*%
      t(centered.measurement)
    gain <- t(solve(S, t(cross.cov)))

    innovation.t <- Y[t, ] - y.pred
    x.filtered <- x.pred + gain %*% innovation.t
    P.filtered <- P.pred - gain %*% S %*% t(gain)
    P.filtered <- (P.filtered + t(P.filtered)) / 2

    predicted.mean[t, ] <- x.pred
    predicted.cov[, , t] <- P.pred
    filtered.mean[t, ] <- x.filtered
    filtered.cov[, , t] <- P.filtered
    innovation[t, ] <- innovation.t
    innovation.cov[, , t] <- S
    loglik.vector[t] <- mvtnorm::dmvnorm(
      innovation.t, sigma = S, log = TRUE
    )

    x.previous <- x.filtered
    P.previous <- P.filtered
  }

  list(
    filtered.mean = filtered.mean,
    filtered.cov = filtered.cov,
    predicted.mean = predicted.mean,
    predicted.cov = predicted.cov,
    innovation = innovation,
    innovation.cov = innovation.cov,
    loglik = sum(loglik.vector),
    loglik.vector = loglik.vector
  )
}
