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

