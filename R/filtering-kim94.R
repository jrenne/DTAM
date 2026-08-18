#' Markov-switching state-space simulation and Kim's filter
#'
#' `simul_SSRS()` simulates a linear Gaussian state-space model whose
#' coefficients depend on a finite-state Markov chain. `Kim94()` applies the
#' collapsing filter of Kim (1994), combining a Kalman update for every
#' transition between regimes with a Hamilton update of regime probabilities.
#'
#' The model is
#' \deqn{y_t = H_{s_t}\beta_t + A_{s_t}Z_t + e_t,\quad
#' e_t \sim N(0,R_{s_t}),}
#' \deqn{\beta_t = \mu_{s_t} + F_{s_t}\beta_{t-1} + G_{s_t}v_t,\quad
#' v_t \sim N(0,Q_{s_t}),}
#' where `s_t` has transition matrix `P` (rows sum to one).
#'
#' @param SSmodel A list containing regime-specific arrays `H`, `A`, `R`,
#'   `mu`, `F`, `G`, and `Q`, plus the regime transition matrix `P`. The third
#'   array dimension indexes regimes.
#' @param s0 Optional initial regime for simulation. If omitted, the most likely
#'   regime under the stationary distribution is used.
#' @param beta0 Initial state mean.
#' @param Z,Z.matrix Matrix of exogenous variables, with one row per date.
#' @param Y.matrix Matrix of observations, with one row per date.
#' @param P0 Initial state covariance matrix.
#' @param xi0 Optional initial regime-probability vector. If omitted, the
#'   stationary distribution of `SSmodel$P` is used.
#'
#' @return `simul_SSRS()` returns simulated regime indicators (`s.simul`),
#'   observations (`y.simul`), and states (`beta.simul`). `Kim94()` returns the
#'   sample log-likelihood (`loglik`), its date-specific contributions
#'   (`loglik.vector`), filtered state means (`all.beta.tt`), filtered state
#'   covariances in vectorized form (`all.P.tt`), and filtered regime
#'   probabilities (`all.xi.t`).
#'
#' @references
#' Kim, C.-J. (1994). Dynamic linear models with Markov-switching.
#' *Journal of Econometrics*, 60, 1--22.
#'
#' Kim, C.-J., and Nelson, C. R. (1999). *State-Space Models with Regime
#' Switching*. MIT Press.
#'
#' @examples
#' # A scalar state with two regimes. Regime 2 has a higher state intercept.
#' model <- list(
#'   H = array(1, c(1, 1, 2)),
#'   A = array(0, c(1, 1, 2)),
#'   R = array(c(0.04, 0.04), c(1, 1, 2)),
#'   mu = array(c(-0.10, 0.10), c(1, 1, 2)),
#'   F = array(c(0.8, 0.8), c(1, 1, 2)),
#'   G = array(1, c(1, 1, 2)),
#'   Q = array(c(0.02, 0.02), c(1, 1, 2)),
#'   P = matrix(c(0.95, 0.05,
#'                0.10, 0.90), 2, 2, byrow = TRUE)
#' )
#' set.seed(123)
#' Z <- matrix(1, 40, 1)
#' sim <- simul_SSRS(model, beta0 = 0, Z = Z)
#' fit <- Kim94(model, sim$y.simul, Z, beta0 = 0, P0 = matrix(0.1))
#' head(fit$all.xi.t)
#' all(abs(rowSums(fit$all.xi.t) - 1) < 1e-10)
#'
#' @name Kim94
NULL

#' @rdname Kim94
#' @export
simul_SSRS <- function(SSmodel, s0 = NaN, beta0, Z) {
  Z <- as.matrix(Z)
  TT <- nrow(Z)
  if (TT < 1L) {
    stop("Z must contain at least one row.")
  }

  H <- SSmodel$H
  A <- SSmodel$A
  R <- SSmodel$R
  mu <- SSmodel$mu
  F <- SSmodel$F
  G <- SSmodel$G
  Q <- SSmodel$Q
  P <- SSmodel$P

  N <- dim(H)[1]
  J <- dim(F)[1]
  L <- dim(G)[2]

  s.simul <- simul_RS(P, TT = TT, ini_state = s0)
  y.simul <- matrix(NA_real_, TT, N)
  beta.simul <- matrix(NA_real_, TT, J)
  beta <- matrix(beta0, J, 1)

  for (t in seq_len(TT)) {
    s <- which.max(s.simul[t, ])
    v <- t(chol(matrix(Q[, , s], L, L))) %*% rnorm(L)
    beta <- mu[, 1, s] + F[, , s] %*% beta + G[, , s] %*% v

    e <- t(chol(matrix(R[, , s], N, N))) %*% rnorm(N)
    y <- matrix(H[, , s], N, J) %*% beta +
      A[, , s] %*% matrix(Z[t, ], ncol = 1) + e

    y.simul[t, ] <- y
    beta.simul[t, ] <- beta
  }

  list(s.simul = s.simul, y.simul = y.simul, beta.simul = beta.simul)
}

#' @rdname Kim94
#' @export
Kim94 <- function(SSmodel, Y.matrix, Z.matrix, beta0, P0, xi0 = NaN) {
  Y.matrix <- as.matrix(Y.matrix)
  Z.matrix <- as.matrix(Z.matrix)
  if (nrow(Y.matrix) != nrow(Z.matrix) || nrow(Y.matrix) < 1L) {
    stop("Y.matrix and Z.matrix must have the same positive number of rows.")
  }

  H <- SSmodel$H
  A <- SSmodel$A
  R <- SSmodel$R
  mu <- SSmodel$mu
  F <- SSmodel$F
  G <- SSmodel$G
  Q <- SSmodel$Q
  P <- as.matrix(SSmodel$P)

  N <- dim(H)[1]
  J <- dim(F)[1]
  M <- dim(H)[3]
  TT <- nrow(Y.matrix)

  if (any(dim(P) != c(M, M)) || any(abs(rowSums(P) - 1) > 1e-10) ||
      any(P < 0)) {
    stop("SSmodel$P must be an M by M transition matrix with rows summing to one.")
  }

  beta.by.regime <- matrix(beta0, J, M)
  P.by.regime <- array(rep(as.matrix(P0), M), c(J, J, M))

  if (is.na(xi0[1])) {
    stationary.system <- rbind(t(P) - diag(M), rep(1, M))
    xi <- as.numeric(qr.solve(stationary.system, c(rep(0, M), 1)))
    xi <- pmax(xi, 0)
    xi <- xi / sum(xi)
  } else {
    xi <- as.numeric(xi0)
    if (length(xi) != M || any(xi < 0) || abs(sum(xi) - 1) > 1e-10) {
      stop("xi0 must be a non-negative probability vector of length M.")
    }
  }

  all.beta.tt <- matrix(NA_real_, TT, J)
  all.P.tt <- matrix(NA_real_, TT, J^2)
  all.xi.t <- matrix(NA_real_, TT, M)
  loglik.vector <- numeric(TT)

  for (t in seq_len(TT)) {
    beta.pair <- array(NA_real_, c(J, M, M))
    P.pair <- array(NA_real_, c(J, J, M, M))
    log.density <- matrix(NA_real_, M, M)

    y.t <- matrix(Y.matrix[t, ], N, 1)
    z.t <- matrix(Z.matrix[t, ], ncol = 1)

    for (i in seq_len(M)) {
      for (j in seq_len(M)) {
        H.j <- matrix(H[, , j], N, J)
        beta.pred <- mu[, 1, j] + F[, , j] %*% beta.by.regime[, i]
        P.pred <- F[, , j] %*% P.by.regime[, , i] %*% t(F[, , j]) +
          G[, , j] %*% Q[, , j] %*% t(G[, , j])
        innovation <- y.t - H.j %*% beta.pred - A[, , j] %*% z.t
        innovation.cov <- H.j %*% P.pred %*% t(H.j) + R[, , j]
        gain <- P.pred %*% t(H.j) %*% solve(innovation.cov)

        beta.pair[, j, i] <- beta.pred + gain %*% innovation
        P.updated <- P.pred - gain %*% H.j %*% P.pred
        P.pair[, , j, i] <- (P.updated + t(P.updated)) / 2
        log.density[i, j] <- mvtnorm::dmvnorm(
          c(innovation), sigma = innovation.cov, log = TRUE
        )
      }
    }

    log.joint <- log.density + log(xi %o% rep(1, M)) + log(P)
    largest <- max(log.joint)
    normalizer <- largest + log(sum(exp(log.joint - largest)))
    posterior.joint <- exp(log.joint - normalizer)
    xi <- colSums(posterior.joint)
    loglik.vector[t] <- normalizer

    next.beta <- matrix(NA_real_, J, M)
    next.P <- array(NA_real_, c(J, J, M))
    for (j in seq_len(M)) {
      if (xi[j] > 0) {
        weights <- posterior.joint[, j] / xi[j]
      } else {
        # This regime is unreachable at date t. Its conditional moments do not
        # enter the aggregate filter; retain a finite placeholder for the next
        # recursion rather than propagating 0 / 0 values.
        weights <- rep(0, M)
        weights[which.max(xi)] <- 1
      }
      next.beta[, j] <- beta.pair[, j, ] %*% weights
      next.P[, , j] <- 0
      for (i in seq_len(M)) {
        difference <- matrix(beta.pair[, j, i] - next.beta[, j], J, 1)
        next.P[, , j] <- next.P[, , j] + weights[i] *
          (P.pair[, , j, i] + difference %*% t(difference))
      }
    }

    beta.filtered <- next.beta %*% xi
    P.filtered <- matrix(0, J, J)
    for (j in seq_len(M)) {
      difference <- matrix(next.beta[, j] - beta.filtered, J, 1)
      P.filtered <- P.filtered + xi[j] *
        (next.P[, , j] + difference %*% t(difference))
    }

    beta.by.regime <- next.beta
    P.by.regime <- next.P
    all.beta.tt[t, ] <- beta.filtered
    all.P.tt[t, ] <- c(P.filtered)
    all.xi.t[t, ] <- xi
  }

  list(
    loglik = sum(loglik.vector),
    loglik.vector = loglik.vector,
    all.beta.tt = all.beta.tt,
    all.P.tt = all.P.tt,
    all.xi.t = all.xi.t
  )
}
