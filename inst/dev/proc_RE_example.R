# Example script for `solve_RE()`.
# This was previously embedded at the end of `R/proc_RE.R`, which made package
# source files harder to maintain. Keeping it here preserves the reference
# calibration without affecting package load behavior.

# Parameter values from Table 1
delta  <- 0.611
kappa  <- 0.064
mu     <- 0.424
rho    <- 0.723
beta   <- 1.525
gamma  <- 0.001
phi_yn <- 0.958
omega  <- 0.877
phi    <- 0.134
phi_1  <- 0.500
phi_2  <- 0.492
phi_3  <- 0.0004
sigma_AS <- 1.249
sigma_IS <- 0.671
sigma_MP <- 2.177
sigma_yn <- 1.380
sigma_pistar <- 0.730

alpha_IS <- 0
alpha_MP <- 0

n <- 5

Gamma0 <- matrix(0, 5, 5)
Gamma0[1, ] <- c(1, -kappa, 0,  kappa, 0)
Gamma0[2, ] <- c(0, 1, phi, 0, 0)
Gamma0[3, ] <- c(0, -(1 - rho) * gamma, 1, (1 - rho) * gamma, (1 - rho) * beta)
Gamma0[4, ] <- c(0, 0, 0, 1, 0)
Gamma0[5, ] <- c(-phi_3, 0, 0, 0, 1)

Gamma1 <- matrix(0, 5, 5)
Gamma1[1, 1] <- 1 - delta
Gamma1[2, 2] <- 1 - mu
Gamma1[3, 3] <- rho
Gamma1[4, 4] <- phi_yn
Gamma1[5, 5] <- phi_2

C <- matrix(c(0, alpha_IS, alpha_MP, 0, 0), 5, 1)

Psi0 <- matrix(0, 5, 5)
Psi0[1, 1] <- delta
Psi0[2, 1] <- phi
Psi0[2, 2] <- mu
Psi0[3, 1] <- (1 - rho) * beta
Psi0[5, 5] <- phi_1

Psi1 <- 0 * Psi0

Pi <- diag(c(sigma_AS, sigma_IS, sigma_MP, sigma_yn, sigma_pistar))

model <- list(
  mu = solve(Gamma0) %*% C,
  Phi = solve(Gamma0) %*% Gamma1,
  Psi0 = solve(Gamma0) %*% Psi0,
  Psi1 = solve(Gamma0) %*% Psi1,
  Sigma12 = solve(Gamma0) %*% Pi
)

res <- solve_RE(model)
