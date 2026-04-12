#' Solve a learning model with filtering
#'
#' Solves a linear learning model with latent states and noisy observations,
#' computes the steady-state filtering objects, and returns the implied joint
#' dynamics of the state vector and its filtered estimate.
#'
#' The model is of the form
#' \deqn{
#' w_t = \mu + \Phi w_{t-1} + \Lambda_0 w_{t|t} + \Lambda_1 w_{t-1|t-1}
#'       + \Sigma^{1/2}\varepsilon_t.
#' }
#'
#' @param model A list containing the model matrices:
#'   `mu`, `Phi`, `Sigma12`, `Omega12`, `A`, `B`, `Lambda0`, and `Lambda1`.
#' @param max.iter Maximum number of iterations used in the fixed-point loop.
#'
#' @return A list containing the original model inputs together with the solved
#'   objects `R`, `P`, `K`, `S`, `mu_ww`, `Phi_ww`, `Sigma_ww`, and
#'   `Sigma12_ww`.
#'
#' @details
#' `Sigma12` and `Omega12` are interpreted as square-root matrices such that
#' `Sigma = Sigma12 %*% t(Sigma12)` and `Omega = Omega12 %*% t(Omega12)`.
#'
#' The function computes the steady-state Kalman gain and then constructs the
#' joint dynamics of the stacked vector
#' \eqn{(w_t^\prime, w_{t|t}^\prime)^\prime}{(w_t', w_{t|t}')'}.
#'
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#'
#' @examples
#' # This function requires a fully specified model list.
#' # See the package vignettes or examples in the source code.
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
#' @export
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
#' Solves an extended learning model in which the state dynamics depend on
#' filtered beliefs and on forward-looking expectations under rational
#' expectations.
#'
#' The model is of the form
#' \deqn{
#' w_t = \mu + \Phi w_{t-1}
#'       + \Lambda_0 w_{t|t} + \Lambda_1 w_{t-1|t-1}
#'       + \Psi_0 E(w_{t+1}|w_t) + \Psi_1 E(w_t|w_{t-1})
#'       + \Gamma_0 E(w_{t+1}|y_t) + \Gamma_1 E(w_t|y_{t-1})
#'       + \Sigma^{1/2}\varepsilon_t.
#' }
#'
#' @param model A list containing the model matrices:
#'   `mu`, `Phi`, `Sigma12`, `Omega12`, `A`, `B`, `Lambda0`, `Lambda1`,
#'   `Psi0`, `Psi1`, `Gamma0`, and `Gamma1`.
#' @param max.iter Maximum number of iterations used in the fixed-point loops.
#'
#' @return A list containing the solved model objects, including `R`, `P`, `K`,
#'   `S`, `mu_ww`, `Phi_ww`, `Sigma_ww`, and `Sigma12_ww`.
#'
#' @details
#' The function first solves for the forward-looking matrix `Phi_1`, then for
#' the feedback matrix `Phi_2`, then computes the steady-state filtering
#' matrices, and finally constructs the joint law of `(w_t', w_{t|t}')'`.
#'
#' A warning is issued if the resulting transition matrix is not stationary.
#'
#' @examples
#' # This function requires a fully specified model list.
#' # See the package vignettes or examples in the source code.
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
