solve_RE <- function(Gamma0, Gamma1, Psi0, Psi1, C, Pi) {

  ## ---- Dimensions and checks ----
  Gamma0 <- as.matrix(Gamma0)
  Gamma1 <- as.matrix(Gamma1)
  Psi0    <- as.matrix(Psi0)
  Psi1    <- as.matrix(Psi1)
  C      <- as.matrix(C)
  Pi     <- as.matrix(Pi)

  n <- nrow(Gamma0)
  if (ncol(Gamma0) != n) stop("Gamma0 must be square.")
  if (!all(dim(Gamma1) == c(n, n))) stop("Gamma1 must be n x n.")
  if (!all(dim(Psi0)   == c(n, n))) stop("Psi0 must be n x n.")
  if (!all(dim(Psi1)   == c(n, n))) stop("Psi1 must be n x n.")
  if (!(nrow(C) == n && ncol(C) == 1)) stop("C must be n x 1.")
  if (nrow(Pi) != n) stop("Pi must have n rows.")

  ## ---- Build the companion pencil ----
  ##   A * s_t = B * s_{t-1}
  ## Characteristic polynomial: Psi*λ^2 - Gamma0*λ + Gamma1 = 0
  A <- rbind(
    cbind(Gamma0, -Psi0),
    cbind(diag(n), 0*diag(n))
  ) + .00000001*diag(2*n)
  B <- rbind(
    cbind(Gamma1, Psi1),
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
  G <- Z2 %*% solve(Z1)

  ## ---- Shock loading and constant term ----
  H <- solve(Gamma0 - Psi0 %*% G) %*% Pi
  d <- solve(Gamma0 - Psi1 - Psi0 %*% (diag(n) + G)) %*% C

  ## ---- Check solution ----
  max_error <- max(abs(Gamma0 %*% G - Gamma1 - Psi0 %*% G %*% G - Psi1 %*% G))

  list(
    G = Re(G),
    H = H,
    d = d,
    max_error = max_error
  )
}

