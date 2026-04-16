## Rational-expectations solver based on a QZ decomposition.
#'
#' Solve a linear rational-expectations model
#'
#' Solves a forward-looking linear system of the form
#' \deqn{
#' w_t = \mu + \Phi w_{t-1} + \Psi_0 E_t(w_{t+1}) + \Psi_1 E_{t-1}(w_t)
#'       + \Sigma^{1/2}\varepsilon_t,
#' }{
#' w_t = mu + Phi w_{t-1} + Psi0 E_t(w_{t+1}) + Psi1 E_{t-1}(w_t) +
#'   Sigma^{1/2} eps_t,
#' }
#' and returns the reduced-form state dynamics implied by the stable
#' rational-expectations solution.
#'
#' @param model A list containing `mu` (state intercept), `Phi`
#'   (lagged-state matrix), `Psi0` and `Psi1` (forward-looking coefficient
#'   matrices), and `Sigma12` (square root of the shock covariance matrix).
#' @return A list with four entries:
#'   `Phi_tilde` (the solved state-transition matrix), `Sigma12_tilde`
#'   (the implied shock-loading matrix), `mu_tilde` (the solved intercept),
#'   and `max_error` (the maximum absolute residual in the matrix fixed-point
#'   equation that characterizes the solution).
#' @details
#' The solution is obtained through a generalized Schur decomposition
#' (QZ decomposition) using [geigen::gqz()]. The function selects the stable
#' generalized eigenvalues and reconstructs the law of motion for `w_t`.
#'
#' `Sigma12` is interpreted as a square-root matrix, so that
#' `Sigma = Sigma12 %*% t(Sigma12)`.
#'
#' This function focuses on the stable rational-expectations solution of the
#' transition equation only. It does not include filtering or measurement
#' equations; those are handled by [solve_learning()] and
#' [solve_Learning_RE()].
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#' @examples
#' # Example adapted from the Bookdown chapter on imperfect information.
#' # It solves a five-dimensional New Keynesian-style rational-expectations
#' # system and returns the stable reduced form.
#' delta <- 0.611
#' kappa <- 0.064
#' mu <- 0.424
#' rho <- 0.723
#' beta <- 1.525
#' gamma <- 0.001
#' phi_yn <- 0.958
#' phi <- 0.134
#' phi_1 <- 0.500
#' phi_2 <- 0.492
#' phi_3 <- 0.0004
#'
#' sigma_AS <- 1.249
#' sigma_IS <- 0.671
#' sigma_MP <- 2.177
#' sigma_yn <- 1.380
#' sigma_pistar <- 0.730
#'
#' alpha_IS <- 0
#' alpha_MP <- 0
#' n <- 5
#'
#' Gamma0 <- matrix(0, 5, 5)
#' Gamma0[1, ] <- c(1, -kappa, 0,  kappa, 0)
#' Gamma0[2, ] <- c(0, 1, phi, 0, 0)
#' Gamma0[3, ] <- c(0, -(1 - rho) * gamma, 1, (1 - rho) * gamma,
#'                  (1 - rho) * beta)
#' Gamma0[4, ] <- c(0, 0, 0, 1, 0)
#' Gamma0[5, ] <- c(-phi_3, 0, 0, 0, 1)
#'
#' Gamma1 <- matrix(0, 5, 5)
#' Gamma1[1, 1] <- 1 - delta
#' Gamma1[2, 2] <- 1 - mu
#' Gamma1[3, 3] <- rho
#' Gamma1[4, 4] <- phi_yn
#' Gamma1[5, 5] <- phi_2
#'
#' C <- matrix(c(0, alpha_IS, alpha_MP, 0, 0), 5, 1)
#'
#' Psi0 <- matrix(0, 5, 5)
#' Psi0[1, 1] <- delta
#' Psi0[2, 1] <- phi
#' Psi0[2, 2] <- mu
#' Psi0[3, 1] <- (1 - rho) * beta
#' Psi0[5, 5] <- phi_1
#'
#' Psi1 <- matrix(0, 5, 5)
#' Pi <- diag(c(sigma_AS, sigma_IS, sigma_MP, sigma_yn, sigma_pistar))
#'
#' model <- list(
#'   mu = solve(Gamma0) %*% C,
#'   Phi = solve(Gamma0) %*% Gamma1,
#'   Psi0 = solve(Gamma0) %*% Psi0,
#'   Psi1 = solve(Gamma0) %*% Psi1,
#'   Sigma12 = solve(Gamma0) %*% Pi
#' )
#'
#' sol <- solve_RE(model)
#' sol$Phi_tilde
#' sol$max_error
#'
#' @export
solve_RE <- function(model) {

  ## ---- Dimensions and checks ----$

  mu      <- as.matrix(model$mu,ncol=1)
  Phi     <- as.matrix(model$Phi)
  Sigma12 <- as.matrix(model$Sigma12)
  Psi0    <- as.matrix(model$Psi0)
  Psi1    <- as.matrix(model$Psi1)

  n <- nrow(Phi)
  if (!all(dim(Phi) == c(n, n))) stop("Phi must be n x n.")
  if (!all(dim(Psi0)   == c(n, n))) stop("Psi0 must be n x n.")
  if (!all(dim(Psi1)   == c(n, n))) stop("Psi1 must be n x n.")
  if (!(nrow(mu) == n && ncol(mu) == 1)) stop("mu must be n x 1.")
  if (nrow(Sigma12) != n) stop("Sigma12 must have n rows.")

  Id_n <- diag(n)

  ## ---- Build the companion pencil ----
  ##   A * s_t = B * s_{t-1}
  ## Characteristic polynomial: Psi*λ^2 - Id_n*λ + Phi = 0
  A <- rbind(
    cbind(Id_n, -Psi0),
    cbind(diag(n), 0*diag(n))
  ) + .00000001*diag(2*n)
  B <- rbind(
    cbind(Phi, Psi1),
    cbind(matrix(0, n, n), diag(n))
  ) + .00000001*diag(2*n)

  qz <- geigen::gqz(A, B, sort="B")

  S <- qz$S
  T <- qz$T
  Q <- qz$Q
  Z <- qz$Z

  # S_ss <- S[1:n,1:n]
  # S_su <- S[1:n,(n+1):(n+n)]
  # S_uu <- S[(n+1):(n+n),(n+1):(n+n)]
  # T_ss <- T[1:n,1:n]
  # T_su <- T[1:n,(n+1):(n+n)]
  # T_uu <- T[(n+1):(n+n),(n+1):(n+n)]
  # Q_ss <- Q[1:n,1:n]
  # Q_su <- Q[1:n,(n+1):(n+n)]
  # Q_us <- Q[(n+1):(n+n),1:n]
  # Q_uu <- Q[(n+1):(n+n),(n+1):(n+n)]
  # Z_ss <- Z[1:n,1:n]
  # Z_su <- Z[1:n,(n+1):(n+n)]
  # Z_us <- Z[(n+1):(n+n),1:n]
  # Z_uu <- Z[(n+1):(n+n),(n+1):(n+n)]
  # print(solve(Z_ss) %*% Z_us)
  # print(solve(S_ss) %*% T_ss)

  # Check QZ decomposition:
  max(abs((A - Q %*% S %*% t(Z))))
  max(abs((B - Q %*% T %*% t(Z))))

  Z_s <- matrix(Z[,1:n],2*n,n)

  ## Partition: Z_s = [ Z1 ; Z2 ]
  Z1 <- matrix(Z_s[1:n,],n,n)
  Z2 <- matrix(Z_s[(n+1):(2*n),],n,n)

  ## ---- Compute G = Z1 %*% solve(Z2) ----
  Phi_tilde <- Z2 %*% solve(Z1)

  ## ---- Shock loading and constant term ----
  Sigma12_tilde <- solve(Id_n - Psi0 %*% Phi_tilde) %*% Sigma12
  mu_tilde <- solve(Id_n - Psi1 - Psi0 %*% (diag(n) + Phi_tilde)) %*% mu

  ## ---- Check solution ----
  max_error <- max(abs(Id_n %*% Phi_tilde - Phi -
                         Psi0 %*% Phi_tilde %*% Phi_tilde - Psi1 %*% Phi_tilde))

  list(
    Phi_tilde = Re(Phi_tilde),
    Sigma12_tilde = Sigma12_tilde,
    mu_tilde = mu_tilde,
    max_error = max_error
  )
}




# Parameter values from Table 1
delta  <- 0.611     # δ
kappa  <- 0.064     # κ
#sigma  <- 3.156     #
mu     <- 0.424     # µ
rho    <- 0.723     # ρ (policy smoothing)
beta   <- 1.525     # β (in policy rule)
gamma  <- 0.001     # γ (response to output gap in policy rule)
phi_yn <- 0.958
omega  <- 0.877     # ω
phi    <- 0.134
phi_1  <- 0.500     # φ₁
phi_2  <- 0.492     # φ₂
phi_3  <- 0.0004 #1 - phi_1 - phi_2  # φ₃ = 1−φ₁−φ₂ ≈ 0.008
sigma_AS  <- 1.249
sigma_IS  <- 0.671
sigma_MP  <- 2.177
sigma_yn  <- 1.380
sigma_pistar <- 0.730

alpha_IS <- 0
alpha_MP <- 0

n <- 5

# Build matrices
Gamma0 <- matrix(0, 5, 5)
Gamma0[1, ] <- c(1, -kappa, 0,  kappa, 0)
Gamma0[2, ] <- c(0,    1,  phi, 0,     0)
Gamma0[3, ] <- c(0, -(1-rho)*gamma, 1, (1-rho)*gamma, (1-rho)*beta)
Gamma0[4, ] <- c(0, 0, 0, 1, 0)
Gamma0[5, ] <- c(-phi_3, 0, 0, 0, 1)

Gamma1 <- matrix(0, 5, 5)
Gamma1[1,1] <- 1 - delta
Gamma1[2,2] <- 1 - mu
Gamma1[3,3] <- rho
Gamma1[4,4] <- phi_yn      # your varphi_{y^n}
Gamma1[5,5] <- phi_2

C <- matrix(c(0, alpha_IS, alpha_MP, 0, 0), 5, 1)

Psi0 <- matrix(0, 5, 5)
Psi0[1,1] <- delta
Psi0[2,1] <- phi
Psi0[2,2] <- mu
Psi0[3,1] <- (1-rho)*beta
Psi0[5,5] <- phi_1

Psi1 <- 0*Psi0

Pi <- diag(c(sigma_AS, sigma_IS, sigma_MP, sigma_yn, sigma_pistar))

## Solve
## Your solver:
model <- list(mu = solve(Gamma0) %*% C,
              Phi = solve(Gamma0) %*% Gamma1,
              Psi0 = solve(Gamma0) %*% Psi0,
              Psi1 = solve(Gamma0) %*% Psi1,
              Sigma12 = solve(Gamma0) %*% Pi)
res <- solve_RE(model)
