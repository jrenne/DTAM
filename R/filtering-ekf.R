as_filter_matrix <- function(x, nrow, ncol, name) {
  x <- as.numeric(x)
  if (length(x) != nrow * ncol || any(!is.finite(x))) {
    stop(name, " must contain ", nrow * ncol, " finite entries.")
  }
  matrix(x, nrow = nrow, ncol = ncol)
}

as_measurement_hessians <- function(x, state.dimension,
                                    measurement.dimension) {
  if (is.list(x)) {
    if (length(x) != measurement.dimension) {
      stop("measurement.hessian must return one Hessian per measurement.")
    }
    hessians <- array(
      NA_real_,
      c(state.dimension, state.dimension, measurement.dimension)
    )
    for (i in seq_len(measurement.dimension)) {
      hessians[, , i] <- as_filter_matrix(
        x[[i]], state.dimension, state.dimension,
        paste0("measurement Hessian ", i)
      )
    }
  } else {
    x <- as.numeric(x)
    expected <- state.dimension^2 * measurement.dimension
    if (length(x) != expected || any(!is.finite(x))) {
      stop(
        "measurement.hessian must return an n by n by m array ",
        "or a list of m n by n matrices."
      )
    }
    hessians <- array(
      x,
      c(state.dimension, state.dimension, measurement.dimension)
    )
  }

  for (i in seq_len(measurement.dimension)) {
    hessians[, , i] <-
      (hessians[, , i] + t(hessians[, , i])) / 2
  }
  hessians
}

#' First- and second-order Extended Kalman filters
#'
#' Applies the first-order (`EKF1`) or second-order (`EKF2`) Extended Kalman
#' filter to a state-space model with an affine transition mean and a smooth,
#' nonlinear measurement function. `EKF` exposes both variants through its
#' `order` argument; `EKF1` and `EKF2` are convenience wrappers.
#'
#' The model is
#' \deqn{x_t = f(x_{t-1},t) + \varepsilon_t,\quad
#' \varepsilon_t \sim N(0,Q),}
#' \deqn{y_t = h(x_t,t) + \eta_t,\quad \eta_t \sim N(0,R).}
#' The transition function `f` must be affine in the state. Its Jacobian is
#' supplied because it propagates the state covariance. EKF1 linearizes `h`.
#' EKF2 additionally uses the Hessian of each measurement to correct the
#' predictive measurement mean and covariance.
#'
#' @param Y Matrix of observations, with one row per date.
#' @param transition Function `transition(x, t)` returning the conditional mean
#'   of the next state. It must be affine in `x`.
#' @param measurement Function `measurement(x, t)` returning the conditional
#'   mean of the observation.
#' @param transition.jacobian Function returning the state-by-state Jacobian of
#'   `transition` at `(x, t)`.
#' @param measurement.jacobian Function returning the measurement-by-state
#'   Jacobian of `measurement` at `(x, t)`.
#' @param Q State-innovation covariance matrix.
#' @param R Measurement-innovation covariance matrix.
#' @param x0 State mean immediately before the first transition.
#' @param P0 Covariance matrix immediately before the first transition.
#' @param order Approximation order: `1` for EKF1 or `2` for EKF2.
#' @param measurement.hessian For EKF2, a function returning either an
#'   `n`-by-`n`-by-`m` array or a list of `m` Hessian matrices, where `n` is the
#'   state dimension and `m` the measurement dimension. For a scalar
#'   measurement, an `n`-by-`n` matrix is sufficient.
#'
#' @return A list containing filtered and predicted state means
#'   (`filtered.mean`, `predicted.mean`), the corresponding covariance arrays
#'   (`filtered.cov`, `predicted.cov`), predictive measurement means
#'   (`predicted.measurement`), innovations and their covariance arrays
#'   (`innovation`, `innovation.cov`), and the total and date-specific Gaussian
#'   quasi-log-likelihood (`loglik`, `loglik.vector`).
#'
#' @details
#' Under the Gaussian predictive approximation, EKF2 adds
#' \deqn{\frac{1}{2}\operatorname{tr}(H_i P)}
#' to the predictive mean of measurement `i`, and
#' \deqn{\frac{1}{2}\operatorname{tr}(H_i P H_j P)}
#' to entry `(i,j)` of the predictive measurement covariance. The
#' state--measurement cross-covariance remains `P A'`, where `A` is the
#' measurement Jacobian.
#'
#' @examples
#' transition <- function(x, t) 0.8 * x
#' transition.jacobian <- function(x, t) matrix(0.8)
#' measurement <- function(x, t) x + 0.1 * x^2
#' measurement.jacobian <- function(x, t) matrix(1 + 0.2 * x)
#' measurement.hessian <- function(x, t) matrix(0.2)
#'
#' set.seed(123)
#' TT <- 30
#' x <- y <- numeric(TT)
#' previous <- 0
#' for (t in seq_len(TT)) {
#'   x[t] <- transition(previous, t) + rnorm(1, sd = 0.2)
#'   y[t] <- measurement(x[t], t) + rnorm(1, sd = 0.1)
#'   previous <- x[t]
#' }
#'
#' fit <- EKF2(
#'   matrix(y), transition, measurement,
#'   transition.jacobian, measurement.jacobian, measurement.hessian,
#'   Q = matrix(0.2^2), R = matrix(0.1^2),
#'   x0 = 0, P0 = matrix(0.2)
#' )
#' head(fit$filtered.mean)
#'
#' @export
EKF <- function(Y, transition, measurement, transition.jacobian,
                measurement.jacobian, Q, R, x0, P0, order = 1L,
                measurement.hessian = NULL) {
  Y <- as.matrix(Y)
  x0 <- as.numeric(x0)
  Q <- as.matrix(Q)
  R <- as.matrix(R)
  P0 <- as.matrix(P0)
  order <- as.integer(order)

  TT <- nrow(Y)
  m <- ncol(Y)
  n <- length(x0)
  if (TT < 1L) {
    stop("Y must contain at least one row.")
  }
  if (anyNA(Y) || any(!is.finite(Y))) {
    stop("EKF does not currently support missing or non-finite observations.")
  }
  if (length(order) != 1L || !order %in% c(1L, 2L)) {
    stop("order must be 1 (EKF1) or 2 (EKF2).")
  }
  if (order == 2L && !is.function(measurement.hessian)) {
    stop("measurement.hessian must be supplied for EKF2.")
  }
  if (any(dim(Q) != c(n, n)) || any(dim(P0) != c(n, n))) {
    stop("Q and P0 must be square matrices matching the state dimension.")
  }
  if (any(dim(R) != c(m, m))) {
    stop("R must be a square matrix matching the observation dimension.")
  }
  if (any(!is.finite(Q)) || any(!is.finite(R)) ||
      any(!is.finite(P0)) || any(!is.finite(x0))) {
    stop("Q, R, x0, and P0 must contain only finite values.")
  }

  Q <- (Q + t(Q)) / 2
  R <- (R + t(R)) / 2
  P0 <- (P0 + t(P0)) / 2

  filtered.mean <- matrix(NA_real_, TT, n)
  predicted.mean <- matrix(NA_real_, TT, n)
  predicted.measurement <- matrix(NA_real_, TT, m)
  innovation <- matrix(NA_real_, TT, m)
  filtered.cov <- array(NA_real_, c(n, n, TT))
  predicted.cov <- array(NA_real_, c(n, n, TT))
  innovation.cov <- array(NA_real_, c(m, m, TT))
  loglik.vector <- numeric(TT)

  x.previous <- x0
  P.previous <- P0

  for (t in seq_len(TT)) {
    F.t <- as_filter_matrix(
      transition.jacobian(x.previous, t), n, n,
      "transition.jacobian"
    )
    x.pred <- as.numeric(transition(x.previous, t))
    if (length(x.pred) != n || any(!is.finite(x.pred))) {
      stop("transition must return a finite vector of state dimension ", n, ".")
    }
    P.pred <- F.t %*% P.previous %*% t(F.t) + Q
    P.pred <- (P.pred + t(P.pred)) / 2

    y.pred <- as.numeric(measurement(x.pred, t))
    if (length(y.pred) != m || any(!is.finite(y.pred))) {
      stop(
        "measurement must return a finite vector of measurement dimension ",
        m, "."
      )
    }
    A.t <- as_filter_matrix(
      measurement.jacobian(x.pred, t), m, n,
      "measurement.jacobian"
    )
    S <- A.t %*% P.pred %*% t(A.t) + R

    if (order == 2L) {
      hessians <- as_measurement_hessians(
        measurement.hessian(x.pred, t), n, m
      )
      for (i in seq_len(m)) {
        H.i <- hessians[, , i]
        y.pred[i] <- y.pred[i] +
          0.5 * sum(diag(H.i %*% P.pred))
        for (j in seq_len(m)) {
          H.j <- hessians[, , j]
          S[i, j] <- S[i, j] + 0.5 *
            sum(diag(H.i %*% P.pred %*% H.j %*% P.pred))
        }
      }
    }

    S <- (S + t(S)) / 2
    cross.cov <- P.pred %*% t(A.t)
    gain <- t(solve(S, t(cross.cov)))

    innovation.t <- Y[t, ] - y.pred
    x.filtered <- c(x.pred + gain %*% innovation.t)
    P.filtered <- P.pred - gain %*% S %*% t(gain)
    P.filtered <- (P.filtered + t(P.filtered)) / 2

    predicted.mean[t, ] <- x.pred
    predicted.cov[, , t] <- P.pred
    predicted.measurement[t, ] <- y.pred
    filtered.mean[t, ] <- x.filtered
    filtered.cov[, , t] <- P.filtered
    innovation[t, ] <- innovation.t
    innovation.cov[, , t] <- S
    loglik.vector[t] <- mvtnorm::dmvnorm(
      innovation.t, mean = rep(0, m), sigma = S, log = TRUE
    )

    x.previous <- x.filtered
    P.previous <- P.filtered
  }

  list(
    filtered.mean = filtered.mean,
    filtered.cov = filtered.cov,
    predicted.mean = predicted.mean,
    predicted.cov = predicted.cov,
    predicted.measurement = predicted.measurement,
    innovation = innovation,
    innovation.cov = innovation.cov,
    loglik = sum(loglik.vector),
    loglik.vector = loglik.vector,
    order = order
  )
}

#' @rdname EKF
#' @export
EKF1 <- function(Y, transition, measurement, transition.jacobian,
                 measurement.jacobian, Q, R, x0, P0) {
  EKF(
    Y, transition, measurement, transition.jacobian,
    measurement.jacobian, Q, R, x0, P0, order = 1L
  )
}

#' @rdname EKF
#' @export
EKF2 <- function(Y, transition, measurement, transition.jacobian,
                 measurement.jacobian, measurement.hessian,
                 Q, R, x0, P0) {
  EKF(
    Y, transition, measurement, transition.jacobian,
    measurement.jacobian, Q, R, x0, P0, order = 2L,
    measurement.hessian = measurement.hessian
  )
}
