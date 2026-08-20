#' Solve a learning model with filtering
#'
#' Solves a linear state-space model with imperfect information, computes the
#' steady-state Kalman filter, and returns a VAR representation for the stacked
#' vector of states and filtered beliefs.
#'
#' The transition and measurement equations are
#' \deqn{
#' w_t = \mu + \Phi w_{t-1} + \Lambda_0 w_{t|t} + \Lambda_1 w_{t-1|t-1}
#'       + \Sigma^{1/2}\varepsilon_t,
#' }{
#' w_t = mu + Phi w_{t-1} + Lambda0 w_{t|t} + Lambda1 w_{t-1|t-1} +
#'   Sigma^{1/2} eps_t,
#' }
#' and
#' \deqn{
#' y_t = B + A w_t + \Omega^{1/2}\eta_t.
#' }{
#' y_t = B + A w_t + Omega^{1/2} eta_t.
#' }
#'
#' @param model A list describing the model. It must contain:
#'   `mu` (state intercept), `Phi` (state-transition matrix), `Sigma12`
#'   (square root of the state-shock covariance matrix), `A` and `B`
#'   (measurement-equation matrices), `Omega12` (square root of the
#'   measurement-error covariance matrix), and `Lambda0`, `Lambda1`
#'   (feedback matrices from filtered beliefs to the transition equation).
#' @param max.iter Maximum number of iterations used to compute the
#'   steady-state filtering covariance matrix.
#' @param tol Convergence tolerance for the covariance recursion. Set to zero
#'   to always run `max.iter` iterations.
#'
#' @return A list containing the original model inputs and the solved objects:
#'   `R`, `P`, `K`, and `S` for the steady-state filter, as well as `mu_ww`,
#'   `Phi_ww`, `Sigma_ww`, and `Sigma12_ww` for the stacked process
#'   \eqn{(w_t^\prime, w_{t|t}^\prime)^\prime}{(w_t', w_{t|t}')'}.
#'   `structural_residuals` reports the direct substitution residuals for the
#'   intercept, transition, and shock equations.
#'
#' @details
#' `Sigma12` and `Omega12` are interpreted as square-root matrices:
#' `Sigma = Sigma12 %*% t(Sigma12)` and
#' `Omega = Omega12 %*% t(Omega12)`.
#'
#' The function iterates on the steady-state Kalman covariance matrix, computes
#' the associated Kalman gain, and then derives the law of motion for the
#' stacked vector formed by the true state and its filtered estimate.
#'
#' The measurement intercept `B` is included in the input specification for
#' consistency with the state-space notation, although the steady-state
#' covariance recursion itself depends only on `A` and `Omega12`.
#'
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#'
#' @examples
#' # Example adapted from the imperfect-information chapter of the companion
#' # Bookdown project. The latent state is two-dimensional and the observation
#' # is noisy consumption growth.
#' set.seed(123)
#'
#' mu_c <- 0
#' phi <- 0.979
#' sigma <- 0.0078
#' varphi_e <- 0.044
#'
#' mu <- matrix(c(0, mu_c), 2, 1)
#' Phi <- matrix(0, 2, 2)
#' Phi[1, 1] <- phi
#' Phi[2, 1] <- 1
#'
#' Sigma12 <- matrix(0, 2, 2)
#' Sigma12[1, 1] <- varphi_e * sigma
#' Sigma12[2, 2] <- sigma
#'
#' A <- matrix(c(0, 1), 1, 2)
#' B <- 0
#' Omega12 <- 1e-7
#'
#' model <- list(
#'   mu = mu,
#'   Phi = Phi,
#'   Sigma12 = Sigma12,
#'   A = A,
#'   B = B,
#'   Omega12 = Omega12,
#'   Lambda0 = 0 * Phi,
#'   Lambda1 = 0 * Phi
#' )
#'
#' model_sol <- solve_learning(model)
#'
#' # Steady-state filter objects:
#' model_sol$K
#' model_sol$P
#'
#' @export
solve_learning <- function(model, max.iter = 200, tol = 1e-10) {

  mu      <- as.matrix(model$mu, ncol = 1)
  Phi     <- as.matrix(model$Phi)
  Sigma12 <- as.matrix(model$Sigma12)
  Omega12 <- as.matrix(model$Omega12)
  A       <- as.matrix(model$A)
  B       <- as.matrix(model$B, ncol = 1)
  Lambda0 <- as.matrix(model$Lambda0)
  Lambda1 <- as.matrix(model$Lambda1)

  Sigma <- Sigma12 %*% t(Sigma12)
  Omega <- Omega12 %*% t(Omega12)

  m <- dim(A)[1]
  n <- dim(A)[2]

  Id_n <- diag(n)

  P0 <- Id_n
  P  <- P0
  converged <- FALSE
  P.change <- Inf

  for (i in seq_len(max.iter)) {
    old_P <- P
    Pstar <- Sigma + Phi %*% P %*% t(Phi)
    S   <- A %*% Pstar %*% t(A) + Omega
    S_1 <- solve(S)
    K   <- Pstar %*% t(A) %*% S_1
    P   <- (Id_n - K %*% A) %*% Pstar
    P.change <- P - old_P
    if (tol > 0 && max(abs(P.change)) <= tol) {
      converged <- TRUE
      break
    }
  }

  R <- solve(Id_n - Lambda0)

  model_sol <- model
  model_sol$R <- R
  model_sol$P <- P
  model_sol$K <- K
  model_sol$S <- S
  model_sol$converged <- converged
  model_sol$iterations <- i
  model_sol$max_change <- max(abs(P.change))

  mu_ww <- rbind(R %*% mu,
                 R %*% mu)

  Theta_1 <- R %*% K %*% A %*% Phi
  Theta_2 <- R %*% ((Id_n - K %*% A) %*% Phi + Lambda1)

  Phi_ww <- rbind(
    cbind(Phi + Lambda0 %*% Theta_1,
          Lambda1 + Lambda0 %*% Theta_2),
    cbind(Theta_1, Theta_2)
  )

  Sigma12_ww <- rbind(
    cbind((Id_n + Lambda0 %*% R %*% K %*% A) %*% Sigma12,
          Lambda0 %*% R %*% K %*% Omega12),
    cbind(R %*% K %*% A %*% Sigma12,
          R %*% K %*% Omega12)
  )

  Sigma_ww <- Sigma12_ww %*% t(Sigma12_ww)

  model_sol$mu_ww      <- mu_ww
  model_sol$Phi_ww     <- Phi_ww
  model_sol$Sigma_ww   <- Sigma_ww
  model_sol$Sigma12_ww <- Sigma12_ww

  top <- seq_len(n)
  bottom <- n + top
  transition_residual <- Phi_ww[top, , drop = FALSE] -
    cbind(Phi, Lambda1) - Lambda0 %*% Phi_ww[bottom, , drop = FALSE]
  shock_residual <- Sigma12_ww[top, , drop = FALSE] -
    cbind(Sigma12, matrix(0, n, ncol(Omega12))) -
    Lambda0 %*% Sigma12_ww[bottom, , drop = FALSE]
  intercept_residual <- mu_ww[top, , drop = FALSE] - mu -
    Lambda0 %*% mu_ww[bottom, , drop = FALSE]
  model_sol$structural_residuals <- c(
    intercept = max(abs(intercept_residual)),
    transition = max(abs(transition_residual)),
    shocks = max(abs(shock_residual))
  )

  return(model_sol)
}

simul_model <- function(model_sol, H) {

  mu_ww      <- model_sol$mu_ww
  Phi_ww     <- model_sol$Phi_ww
  Sigma12_ww <- t(chol(model_sol$Sigma_ww))

  n_ww <- length(mu_ww)
  n    <- n_ww / 2

  x <- solve(diag(n_ww) - Phi_ww) %*% mu_ww
  X <- matrix(x, ncol = 1)

  for (t in 2:H) {
    x <- mu_ww + Phi_ww %*% x + Sigma12_ww %*% rnorm(dim(Sigma12_ww)[2])
    X <- cbind(X, x)
  }

  w_t  <- X[1:n, ]
  w_tt <- X[(n + 1):(2 * n), ]

  return(list(
    w_t  = w_t,
    w_tt = w_tt
  ))
}

#' Solve a learning model with rational-expectations feedback
#'
#' Solves an imperfect-information linear model in which state dynamics depend
#' both on filtered beliefs and on forward-looking rational expectations, and
#' returns the resulting stacked law of motion for states and filtered beliefs.
#'
#' The transition equation is
#' \deqn{
#' w_t = \mu + \Phi w_{t-1}
#'       + \Lambda_0 w_{t|t} + \Lambda_1 w_{t-1|t-1}
#'       + \Psi_0 E(w_{t+1}|w_t) + \Psi_1 E(w_t|w_{t-1})
#'       + \Gamma_0 E(w_{t+1}|y_t) + \Gamma_1 E(w_t|y_{t-1})
#'       + \Sigma^{1/2}\varepsilon_t,
#' }{
#' w_t = mu + Phi w_{t-1} + Lambda0 w_{t|t} + Lambda1 w_{t-1|t-1} +
#'   Psi0 E(w_{t+1}|w_t) + Psi1 E(w_t|w_{t-1}) +
#'   Gamma0 E(w_{t+1}|y_t) + Gamma1 E(w_t|y_{t-1}) +
#'   Sigma^{1/2} eps_t,
#' }
#'
#' The observation equation is
#' \deqn{
#' y_t = B + A w_t + \Omega^{1/2}\eta_t.
#' }{
#' y_t = B + A w_t + Omega^{1/2} eta_t.
#' }
#'
#' @param model A list describing the model. It must contain:
#'   `mu`, `Phi`, `Sigma12`, `Omega12`, `A`, `B`, `Lambda0`, `Lambda1`,
#'   `Psi0`, `Psi1`, `Gamma0`, and `Gamma1`. Here `Lambda0` and `Lambda1`
#'   capture feedback from filtered beliefs, while `Psi0`, `Psi1`, `Gamma0`,
#'   and `Gamma1` govern the effects of forward-looking expectations formed
#'   under rational expectations.
#' @param max.iter Maximum number of iterations used in the fixed-point and
#'   filtering loops.
#' @param tol Convergence tolerance for each fixed-point loop. Set to zero to
#'   always run `max.iter` iterations.
#' @param damping Number in `(0, 1]` used to damp updates of the coupled
#'   transition/filter fixed point. Smaller values are slower but can improve
#'   convergence in models with strong expectational feedback.
#'
#' @return A list containing the original model inputs together with the solved
#'   objects `R`, `P`, `K`, `S`, `mu_ww`, `Phi_ww`, `Sigma_ww`, and
#'   `Sigma12_ww`. The element `structural_residuals` reports the maximum
#'   absolute residual obtained by substituting the solution back into the
#'   original intercept, transition, and shock equations.
#'
#' @details
#' `Sigma12` and `Omega12` are interpreted as square-root matrices:
#' `Sigma = Sigma12 %*% t(Sigma12)` and
#' `Omega = Omega12 %*% t(Omega12)`.
#'
#' The two transition blocks and the steady-state Kalman filter are mutually
#' dependent when contemporaneous filtered expectations enter the structural
#' equations. The function therefore solves them as a coupled damped fixed
#' point and then builds the joint dynamics of the stacked vector
#' \eqn{(w_t^\prime, w_{t|t}^\prime)^\prime}{(w_t', w_{t|t}')'}.
#'
#' The returned list contains the original model inputs together with the
#' solved objects `R`, `P`, `K`, `S`, `mu_ww`, `Phi_ww`, `Sigma_ww`, and
#' `Sigma12_ww`.
#'
#' Warnings are issued if the coupled iteration or filter does not converge,
#' if direct substitution produces non-negligible structural residuals, or if
#' the resulting stacked transition matrix is non-stationary. Convergence of
#' the numerical iteration alone is not treated as an equilibrium check.
#'
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#'
#' @examples
#' model <- list(
#'   mu = matrix(0), Phi = matrix(0.7),
#'   Sigma12 = matrix(0.1), Omega12 = matrix(0.2),
#'   A = matrix(1), B = matrix(0),
#'   Lambda0 = matrix(0.02), Lambda1 = matrix(0.01),
#'   Psi0 = matrix(0.05), Psi1 = matrix(0.02),
#'   Gamma0 = matrix(0.02), Gamma1 = matrix(0.01)
#' )
#'
#' model_sol <- solve_Learning_RE(model)
#' model_sol$K
#' model_sol$Phi_ww
#'
#' @importFrom MASS ginv
#' @export
solve_Learning_RE <- function(model, max.iter = 1000, tol = 1e-10,
                              damping = 0.15) {

  mu      <- as.matrix(model$mu, ncol = 1)
  Phi     <- as.matrix(model$Phi)
  Sigma12 <- as.matrix(model$Sigma12)
  Omega12 <- as.matrix(model$Omega12)
  A       <- as.matrix(model$A)
  B       <- as.matrix(model$B, ncol = 1)
  Lambda0 <- as.matrix(model$Lambda0)
  Lambda1 <- as.matrix(model$Lambda1)
  Psi0    <- as.matrix(model$Psi0)
  Psi1    <- as.matrix(model$Psi1)
  Gamma0  <- as.matrix(model$Gamma0)
  Gamma1  <- as.matrix(model$Gamma1)

  Omega <- Omega12 %*% t(Omega12)

  n <- dim(A)[2]
  Id_n <- diag(n)

  if (length(max.iter) != 1 || max.iter < 1) {
    stop("max.iter must be a positive integer.")
  }
  if (length(tol) != 1 || !is.finite(tol) || tol < 0) {
    stop("tol must be a non-negative number.")
  }
  if (length(damping) != 1 || !is.finite(damping) ||
      damping <= 0 || damping > 1) {
    stop("damping must be a number in (0, 1].")
  }

  # The full-information solution and the former sequential feedback update
  # provide an effective warm start. The equilibrium is then solved from the
  # coupled transition/filter equations below.
  Phi_1 <- Phi
  for (iter.init1 in seq_len(max.iter)) {
    old.Phi_1 <- Phi_1
    Phi_1 <- Phi +
      Psi1 %*% Phi_1 +
      Psi0 %*% Phi_1 %*% Phi_1
    if (tol > 0 && max(abs(Phi_1 - old.Phi_1)) <= tol) {
      break
    }
  }

  Phi_2 <- matrix(0, n, n)
  for (iter.init2 in seq_len(max.iter)) {
    old.Phi_2 <- Phi_2
    M_init <- solve(Id_n - Psi0 %*% Phi_1)
    Lambda0_init <- M_init %*%
      (Psi0 %*% Phi_2 + Lambda0 + Gamma0 %*% (Phi_1 + Phi_2))
    Lambda1_init <- M_init %*%
      (Psi1 %*% Phi_2 + Lambda1 + Gamma1 %*% (Phi_1 + Phi_2))
    R_init <- solve(Id_n - Lambda0_init)
    Phi_2 <- (R_init - Id_n) %*% Phi_1 + R_init %*% Lambda1_init
    if (tol > 0 && max(abs(Phi_2 - old.Phi_2)) <= tol) {
      break
    }
  }

  fixed_point_map <- function(Phi_1, Phi_2, P) {
    M_tilde <- solve(Id_n - Psi0 %*% Phi_1)
    F_tilde <- M_tilde %*% (Phi + Psi1 %*% Phi_1)
    Lambda0_tilde <- M_tilde %*%
      (Psi0 %*% Phi_2 + Lambda0 + Gamma0 %*% (Phi_1 + Phi_2))
    Lambda1_tilde <- M_tilde %*%
      (Psi1 %*% Phi_2 + Lambda1 + Gamma1 %*% (Phi_1 + Phi_2))
    Sigma12_tilde <- M_tilde %*% Sigma12
    Sigma_tilde <- Sigma12_tilde %*% t(Sigma12_tilde)
    R <- solve(Id_n - Lambda0_tilde)

    converged.filter <- FALSE
    P.change <- Inf
    for (iter.filter in seq_len(max.iter)) {
      old_P <- P
      Pstar <- Sigma_tilde + F_tilde %*% P %*% t(F_tilde)
      S <- A %*% Pstar %*% t(A) + Omega
      K <- Pstar %*% t(A) %*% MASS::ginv(S)
      P <- (Id_n - K %*% A) %*% Pstar
      P <- (P + t(P)) / 2
      P.change <- P - old_P
      if (tol > 0 && max(abs(P.change)) <= tol) {
        converged.filter <- TRUE
        break
      }
    }
    if (tol == 0 && max(abs(P.change)) <= 1e-8) {
      converged.filter <- TRUE
    }

    Theta_1 <- R %*% K %*% A %*% F_tilde
    Theta_2 <- R %*% ((Id_n - K %*% A) %*% F_tilde + Lambda1_tilde)

    list(
      Phi_1 = F_tilde + Lambda0_tilde %*% Theta_1,
      Phi_2 = Lambda1_tilde + Lambda0_tilde %*% Theta_2,
      M_tilde = M_tilde,
      F_tilde = F_tilde,
      Lambda0_tilde = Lambda0_tilde,
      Lambda1_tilde = Lambda1_tilde,
      Sigma12_tilde = Sigma12_tilde,
      Sigma_tilde = Sigma_tilde,
      R = R,
      P = P,
      Pstar = Pstar,
      K = K,
      S = S,
      Theta_1 = Theta_1,
      Theta_2 = Theta_2,
      converged.filter = converged.filter,
      iter.filter = iter.filter,
      filter.change = max(abs(P.change))
    )
  }

  P <- Id_n
  converged.joint <- FALSE
  joint.change <- Inf
  for (iter.joint in seq_len(max.iter)) {
    step <- fixed_point_map(Phi_1, Phi_2, P)
    joint.change <- max(abs(c(step$Phi_1 - Phi_1,
                              step$Phi_2 - Phi_2)))
    P <- step$P
    if (tol > 0 && joint.change <= tol) {
      Phi_1 <- step$Phi_1
      Phi_2 <- step$Phi_2
      converged.joint <- TRUE
      break
    }
    Phi_1 <- (1 - damping) * Phi_1 + damping * step$Phi_1
    Phi_2 <- (1 - damping) * Phi_2 + damping * step$Phi_2
  }

  # Recompute every object at the final transition matrices.
  step <- fixed_point_map(Phi_1, Phi_2, P)
  joint.change <- max(abs(c(step$Phi_1 - Phi_1,
                            step$Phi_2 - Phi_2)))
  if ((tol > 0 && joint.change <= tol) ||
      (tol == 0 && joint.change <= 1e-8)) {
    converged.joint <- TRUE
  }

  M_tilde <- step$M_tilde
  F_tilde <- step$F_tilde
  Lambda0_tilde <- step$Lambda0_tilde
  Lambda1_tilde <- step$Lambda1_tilde
  Sigma12_tilde <- step$Sigma12_tilde
  Sigma_tilde <- step$Sigma_tilde
  R <- step$R
  P <- step$P
  K <- step$K
  S <- step$S
  Theta_1 <- step$Theta_1
  Theta_2 <- step$Theta_2

  D <- Psi0 + Psi1 + Gamma0 + Gamma1
  mu_tilde <- solve(
    (Id_n - Psi0 %*% Phi_1) %*% (Id_n - Lambda0_tilde) - D,
    mu
  )

  H11 <- (Id_n + Lambda0_tilde %*% R %*% K %*% A) %*%
    Sigma12_tilde
  H12 <- Lambda0_tilde %*% R %*% K %*% Omega12
  H21 <- R %*% K %*% A %*% Sigma12_tilde
  H22 <- R %*% K %*% Omega12

  model_sol <- model
  model_sol$R <- R
  model_sol$P <- P
  model_sol$K <- K
  model_sol$S <- S
  model_sol$M_tilde <- M_tilde
  model_sol$F_tilde <- F_tilde
  model_sol$Lambda0_tilde <- Lambda0_tilde
  model_sol$Lambda1_tilde <- Lambda1_tilde
  model_sol$Sigma12_tilde <- Sigma12_tilde
  model_sol$Sigma_tilde <- Sigma_tilde
  model_sol$Phi_1 <- Phi_1
  model_sol$Phi_2 <- Phi_2

  mu_ww <- rbind(mu_tilde, mu_tilde)

  Phi_ww <- rbind(
    cbind(Phi_1, Phi_2),
    cbind(Theta_1, Theta_2)
  )

  Sigma12_ww <- rbind(
    cbind(H11, H12),
    cbind(H21, H22)
  )

  Sigma_ww <- Sigma12_ww %*% t(Sigma12_ww)

  top <- seq_len(n)
  bottom <- n + top
  Jw <- cbind(Id_n, matrix(0, n, n))
  Jbelief <- cbind(matrix(0, n, n), Id_n)
  private_transition <- Phi_1 + Phi_2
  transition_rhs <- Phi %*% Jw +
    Lambda0 %*% Phi_ww[bottom, , drop = FALSE] +
    Lambda1 %*% Jbelief +
    Psi0 %*% Phi_ww[top, , drop = FALSE] %*% Phi_ww +
    Psi1 %*% Phi_ww[top, , drop = FALSE] +
    Gamma0 %*% private_transition %*% Phi_ww[bottom, , drop = FALSE] +
    Gamma1 %*% private_transition %*% Jbelief
  transition_residual <- Phi_ww[top, , drop = FALSE] - transition_rhs

  structural_shocks <- cbind(
    Sigma12,
    matrix(0, n, ncol(Omega12))
  )
  shock_rhs <- Lambda0 %*% Sigma12_ww[bottom, , drop = FALSE] +
    Psi0 %*% Phi_ww[top, , drop = FALSE] %*% Sigma12_ww +
    Gamma0 %*% private_transition %*%
      Sigma12_ww[bottom, , drop = FALSE] + structural_shocks
  shock_residual <- Sigma12_ww[top, , drop = FALSE] - shock_rhs

  intercept_rhs <- mu + (Lambda0 + D) %*% mu_tilde +
    (Psi0 + Gamma0) %*% private_transition %*% mu_tilde
  intercept_residual <- mu_tilde - intercept_rhs
  structural_residuals <- c(
    intercept = max(abs(intercept_residual)),
    transition = max(abs(transition_residual)),
    shocks = max(abs(shock_residual))
  )
  check.tol <- if (tol > 0) max(100 * tol, 1e-8) else 1e-8
  converged.structural <- max(structural_residuals) <= check.tol

  model_sol$converged <- c(
    joint = converged.joint,
    filter = step$converged.filter,
    structural = converged.structural
  )
  model_sol$iterations <- c(
    joint = iter.joint,
    filter = step$iter.filter
  )
  model_sol$max_change <- c(
    joint = joint.change,
    filter = step$filter.change,
    structural = max(structural_residuals)
  )
  model_sol$structural_residuals <- structural_residuals

  if (!converged.joint) {
    warning("The coupled transition/filter fixed point did not converge.")
  }
  if (!step$converged.filter) {
    warning("The steady-state filtering covariance did not converge.")
  }
  if (!converged.structural) {
    warning(
      "The reported solution does not satisfy the original structural ",
      "equations within tolerance (maximum residual = ",
      format(max(structural_residuals), digits = 4), ")."
    )
  }

  res_eig <- eigen(Phi_ww)
  if (sum(abs(res_eig$values) >= 1) > 0) {
    warning(
      "The resulting dynamics of w_t is not stationary ",
      "(eigenvalues of Phi_tilde exceeding 1 in modulus)."
    )
  }

  model_sol$mu_ww      <- mu_ww
  model_sol$Phi_ww     <- Phi_ww
  model_sol$Sigma_ww   <- Sigma_ww
  model_sol$Sigma12_ww <- Sigma12_ww

  return(model_sol)
}
