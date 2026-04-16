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
#'
#' @return A list containing the original model inputs and the solved objects:
#'   `R`, `P`, `K`, and `S` for the steady-state filter, as well as `mu_ww`,
#'   `Phi_ww`, `Sigma_ww`, and `Sigma12_ww` for the stacked process
#'   \eqn{(w_t^\prime, w_{t|t}^\prime)^\prime}{(w_t', w_{t|t}')'}.
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
solve_learning <- function(model, max.iter = 200) {

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
  old_P <- P

  for (i in 1:max.iter) {
    Pstar <- Sigma + Phi %*% P %*% t(Phi)
    S   <- A %*% Pstar %*% t(A) + Omega
    S_1 <- solve(S)
    K   <- Pstar %*% t(A) %*% S_1
    P   <- (Id_n - K %*% A) %*% Pstar
    P.change <- P - old_P
    old_P <- P
  }

  R <- solve(Id_n - Lambda0)

  model_sol <- model
  model_sol$R <- R
  model_sol$P <- P
  model_sol$K <- K
  model_sol$S <- S

  mu_ww <- rbind(R %*% mu,
                 R %*% mu)

  Phi_ww <- rbind(
    cbind(Phi,
          Lambda1 + Lambda0 %*% R %*% (Phi + Lambda1)),
    cbind(R %*% K %*% A %*% Phi,
          R %*% ((Id_n - K %*% A) %*% Phi + Lambda1))
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

  return(model_sol)
}

#' Simulate the solved learning model
#'
#' Simulates the joint process of the state vector and its filtered estimate
#' from the solved representation returned by [solve_learning()] or
#' [solve_Learning_RE()].
#'
#' @param model_sol A solved model object containing at least `mu_ww`,
#'   `Phi_ww`, and `Sigma_ww`.
#' @param H Simulation horizon.
#'
#' @return A list with two matrices:
#'   \itemize{
#'   \item `w_t`: simulated state vector,
#'   \item `w_tt`: simulated filtered state vector.
#'   }
#'
#' @examples
#' # Requires a solved model object:
#' # sim <- simul_model(model_sol, H = 200)
#'
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
#' @param max.iter Maximum number of iterations used in the successive
#'   fixed-point loops for the rational-expectations and filtering objects.
#'
#' @return A list containing the original model inputs together with the solved
#'   objects `R`, `P`, `K`, `S`, `mu_ww`, `Phi_ww`, `Sigma_ww`, and
#'   `Sigma12_ww`.
#'
#' @details
#' `Sigma12` and `Omega12` are interpreted as square-root matrices:
#' `Sigma = Sigma12 %*% t(Sigma12)` and
#' `Omega = Omega12 %*% t(Omega12)`.
#'
#' The function first solves for the forward-looking matrix `Phi_1`, then for
#' the feedback matrix `Phi_2`, next computes the steady-state Kalman filter,
#' and finally builds the joint dynamics of the stacked vector
#' \eqn{(w_t^\prime, w_{t|t}^\prime)^\prime}{(w_t', w_{t|t}')'}.
#'
#' The returned list contains the original model inputs together with the
#' solved objects `R`, `P`, `K`, `S`, `mu_ww`, `Phi_ww`, `Sigma_ww`, and
#' `Sigma12_ww`.
#'
#' A warning is issued if the eigenvalues of the resulting stacked transition
#' matrix indicate non-stationary dynamics.
#'
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#'
#' @examples
#' # Example adapted from the imperfect-information chapter of the companion
#' # Bookdown project. It solves a five-dimensional New Keynesian model with
#' # rational expectations and noisy signals on a subset of the states.
#' delta <- 0.611
#' kappa <- 0.064
#' rho <- 0.723
#' beta <- 1.525
#' gamma <- 0.001
#' phi_yn <- 0.958
#' phi_1 <- 0.500
#' phi_2 <- 0.492
#' phi_3 <- 0.0004
#'
#' sigma_AS <- 1.249
#' sigma_IS <- 0.671
#' sigma_MP <- 2.177
#' sigma_rstar <- 1.380
#' sigma_pistar <- 0.730
#' sigma_nu <- 0.5
#'
#' Theta <- matrix(c(
#'   1, -kappa, 0, kappa, 0,
#'   0, 1, 0.134, 0, 0,
#'   0, -(1 - rho) * gamma, 1, (1 - rho) * gamma, (1 - rho) * beta,
#'   0, 0, 0, 1, 0,
#'   -phi_3, 0, 0, 0, 1
#' ), nrow = 5, byrow = TRUE)
#'
#' Theta_inv <- solve(Theta)
#' mu <- matrix(0, nrow = 5, ncol = 1)
#'
#' Phi_raw <- matrix(0, 5, 5)
#' Phi_raw[1, 1] <- 1 - delta
#' Phi_raw[2, 2] <- 1 - 0.424
#' Phi_raw[3, 3] <- rho
#' Phi_raw[4, 4] <- phi_yn
#' Phi_raw[5, 5] <- phi_2
#' Phi <- Theta_inv %*% Phi_raw
#'
#' Gamma0_raw <- matrix(0, 5, 5)
#' Gamma0_raw[1, ] <- c(delta, -kappa, 0, kappa, 0)
#' Gamma0_raw[2, ] <- c(0.134, 0.424, 0, 0, 0)
#' Gamma0_raw[3, ] <- c((1 - rho) * beta, 0, 0, 0, 0)
#' Gamma0_raw[5, 5] <- phi_1
#' Gamma0 <- Theta_inv %*% Gamma0_raw
#'
#' Psi0 <- matrix(0, 5, 5)
#' Lambda0 <- matrix(0, 5, 5)
#' Lambda1 <- matrix(0, 5, 5)
#' Psi1 <- matrix(0, 5, 5)
#' Gamma1 <- matrix(0, 5, 5)
#'
#' Sigma_half_raw <- diag(c(
#'   sigma_AS,
#'   sigma_IS,
#'   sigma_MP,
#'   sigma_rstar,
#'   sigma_pistar
#' ))
#' Sigma12 <- Theta_inv %*% Sigma_half_raw
#'
#' A <- matrix(c(
#'   1, 0, 0, 0, 0,
#'   0, 1, 0, 0, 0,
#'   0, 0, 1, 0, 0
#' ), nrow = 3, byrow = TRUE)
#' B <- matrix(0, nrow = 3, ncol = 1)
#' Omega12 <- matrix(c(0, sigma_nu, 0), nrow = 3, ncol = 1)
#'
#' model <- list(
#'   mu = mu,
#'   Phi = Phi,
#'   Gamma0 = Gamma0,
#'   Psi0 = Psi0,
#'   Sigma12 = Sigma12,
#'   Lambda0 = Lambda0,
#'   Lambda1 = Lambda1,
#'   Psi1 = Psi1,
#'   Gamma1 = Gamma1,
#'   Omega12 = Omega12,
#'   A = A,
#'   B = B
#' )
#'
#' model_sol <- solve_Learning_RE(model)
#' model_sol$K
#' model_sol$Phi_ww
#'
#' @importFrom MASS ginv
#' @export
solve_Learning_RE <- function(model, max.iter = 200) {

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

  Sigma <- Sigma12 %*% t(Sigma12)
  Omega <- Omega12 %*% t(Omega12)

  m <- dim(A)[1]
  n <- dim(A)[2]
  Id_n <- diag(n)

  Phi_1 <- Phi
  for (i in 1:max.iter) {
    Phi_1 <- Phi +
      Psi1 %*% Phi_1 +
      Psi0 %*% Phi_1 %*% Phi_1
  }

  Id_Psi0_Phi1_1 <- solve(Id_n - Psi0 %*% Phi_1)

  Phi_2 <- matrix(0, n, n)
  for (i in 1:max.iter) {

    Lambda0_tilde <- Id_Psi0_Phi1_1 %*%
      (Psi0 %*% Phi_2 + Lambda0 + Gamma0 %*% (Phi_1 + Phi_2))

    Lambda1_tilde <- Id_Psi0_Phi1_1 %*%
      (Psi1 %*% Phi_2 + Lambda1 + Gamma1 %*% (Phi_1 + Phi_2))

    Phi_2 <- (solve(Id_n - Lambda0_tilde) - Id_n) %*% Phi_1 +
      solve(Id_n - Lambda0_tilde) %*% Lambda1_tilde
  }

  R <- solve(Id_n - Lambda0_tilde)

  mu_tilde <- R %*% Id_Psi0_Phi1_1 %*%
    (Id_n + Psi0 + Psi1 + Gamma0 + Gamma1) %*% mu

  Sigma_tilde <- Id_Psi0_Phi1_1 %*% Sigma %*% t(Id_Psi0_Phi1_1)

  P0 <- Id_n
  P  <- P0
  old_P <- P

  for (i in 1:max.iter) {

    Pstar <- Sigma_tilde + Phi_1 %*% P %*% t(Phi_1)

    S   <- A %*% Pstar %*% t(A) + Omega
    S_1 <- MASS::ginv(S)

    K <- Pstar %*% t(A) %*% S_1

    P <- (Id_n - K %*% A) %*% Pstar

    P.change <- P - old_P
    old_P <- P
  }

  Theta_1 <- R %*% K %*% A %*% Phi_1
  Theta_2 <- R %*% ((Id_n - K %*% A) %*% Phi_1 + Lambda1_tilde)

  H11 <- (Id_n + Lambda0_tilde %*% R %*% K %*% A) %*% Sigma12
  H12 <- Lambda0_tilde %*% R %*% K %*% Omega12
  H21 <- R %*% K %*% A %*% Sigma12
  H22 <- R %*% K %*% Omega12

  model_sol <- model
  model_sol$R <- R
  model_sol$P <- P
  model_sol$K <- K
  model_sol$S <- S

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
