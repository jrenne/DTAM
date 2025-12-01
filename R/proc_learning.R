


solve_learning <- function(model,max.iter=200){
  #
  # The model is:
  # w_t = mu + Phi·w_{t+1} + Lambda0·w_{t|t} + Lambda1·w_{t-1|t-1}
  #         + Sigma^{1/2}·epsilon_t
  #

  mu      <- as.matrix(model$mu,ncol=1)
  Phi     <- as.matrix(model$Phi)
  Sigma12 <- as.matrix(model$Sigma12)
  Omega12 <- as.matrix(model$Omega12)
  A       <- as.matrix(model$A)
  B       <- as.matrix(model$B,ncol=1)
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

  for (i in 1:max.iter){
    Pstar <- Sigma + Phi %*% P %*% t(Phi)
    S   <- A %*% Pstar %*% t(A) + Omega
    S_1 <- solve(S)
    K   <- Pstar %*% t(A) %*% S_1
    P  <- (Id_n - Pstar %*% t(A) %*%
             solve(A %*% Pstar %*% t(A)) %*% A) %*% Pstar
    P.change <- P - old_P
    old_P <- P
  }

  R <- solve(Id_n - Lambda0)

  model_sol <- model
  model_sol$R <- R
  model_sol$P <- P
  model_sol$K <- K
  model_sol$S <- S

  # Dynamics of (w_t',w_{t|t}')':

  mu_ww <- rbind(R %*% mu,
                 R %*% mu)
  Phi_ww <- rbind(cbind(Phi,
                        Lambda1 + Lambda0 %*% R %*% (Phi + Lambda1)),
                  cbind(R %*% K %*% A %*% Phi,
                        R %*% ((Id_n - K %*% A) %*% Phi + Lambda1)))

  Sigma12_ww <- rbind(cbind((Id_n + Lambda0 %*% R %*% K %*% A) %*% Sigma12,
                            Lambda0 %*% R %*% K %*% Omega12),
                      cbind(R %*% K %*% A %*% Sigma12,
                            R %*% K %*% Omega12))

  Sigma_ww <- Sigma12_ww %*% t(Sigma12_ww)

  model_sol$mu_ww    <- mu_ww
  model_sol$Phi_ww   <- Phi_ww
  model_sol$Sigma_ww <- Sigma_ww
  model_sol$Sigma12_ww <- Sigma12_ww

  return(model_sol)
}


simul_model <- function(model_sol,H){

  mu_ww      <- model_sol$mu_ww
  Phi_ww     <- model_sol$Phi_ww
  Sigma12_ww <- t(chol(model_sol$Sigma_ww))

  n_ww <- length(mu_ww)
  n   <- n_ww/2

  x <- solve(diag(n_ww) - Phi_ww) %*% mu_ww
  X <- matrix(x,ncol=1)

  for(t in 2:H){
    x <- mu_ww + Phi_ww %*% x + Sigma12_ww %*% rnorm(dim(Sigma12_ww)[2])
    X <- cbind(X,x)
  }
  w_t  <- X[1:n,]
  w_tt <- X[(n+1):(2*n),]

  return(list(w_t  = w_t,
              w_tt = w_tt))
}


solve_Learning_RE <- function(model, max.iter = 200){

  #
  # The model is:
  # w_t = mu + Phi·w_{t+1}
  #          + Lambda0·w_{t|t}         + Lambda1·w_{t-1|t-1}
  #          + Psi0·E(w_{t+1}|w_{t})   + Psi1·E(w_{t}|w_{t-1})
  #          + Gamma0·E(w_{t+1}|y_{t}) + Gamma1·E(w_{t}|y_{t-1})
  #          + Sigma^{1/2}·epsilon_t
  #

  # --- Extract model inputs and ensure matrix format --------------------------

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

  # Construct covariance matrices
  Sigma <- Sigma12 %*% t(Sigma12)
  Omega <- Omega12 %*% t(Omega12)

  # Dimensions
  m <- dim(A)[1]
  n <- dim(A)[2]
  Id_n <- diag(n)

  # ----------------------------------------------------------------------------
  # 1. Solve for Phi_1 (forward-looking expectations)
  # ----------------------------------------------------------------------------

  Phi_1 <- Phi
  for (i in 1:max.iter){
    Phi_1 <- Phi +
      Psi1 %*% Phi_1 +
      Psi0 %*% Phi_1 %*% Phi_1
  }

  # ----------------------------------------------------------------------------
  # 2. Solve recursively for Phi_2 (feedback terms)
  # ----------------------------------------------------------------------------

  Id_Psi0_Phi1_1 <- solve(Id_n - Psi0 %*% Phi_1)

  Phi_2 <- matrix(0, n, n)
  for (i in 1:max.iter){

    # Time-tilded Lambda matrices
    Lambda0_tilde <- Id_Psi0_Phi1_1 %*%
      (Psi0 %*% Phi_2 + Lambda0 + Gamma0 %*% (Phi_1 + Phi_2))

    Lambda1_tilde <- Id_Psi0_Phi1_1 %*%
      (Psi1 %*% Phi_2 + Lambda1 + Gamma1 %*% (Phi_1 + Phi_2))

    # Update Phi_2 using fixed-point iteration
    Phi_2 <- (solve(Id_n - Lambda0_tilde) - Id_n) %*% Phi_1 +
      solve(Id_n - Lambda0_tilde) %*% Lambda1_tilde
  }

  # Inverse used repeatedly
  R <- solve(Id_n - Lambda0_tilde)

  # Compute mu_tilde (constant term)
  mu_tilde <- R %*% Id_Psi0_Phi1_1 %*%
    (Id_n + Psi0 + Psi1 + Gamma0 + Gamma1) %*% mu

  # ----------------------------------------------------------------------------
  # 3. Solve for the Kalman gain K (filtering step)
  # ----------------------------------------------------------------------------

  # “Tilted” Sigma matrix
  Sigma_tilde <- Id_Psi0_Phi1_1 %*% Sigma %*% t(Id_Psi0_Phi1_1)

  # Initialize P
  P0 <- Id_n
  P  <- P0
  old_P <- P

  for (i in 1:max.iter){

    # Prediction variance
    Pstar <- Sigma_tilde + Phi_1 %*% P %*% t(Phi_1)

    # Observation variance
    S   <- A %*% Pstar %*% t(A) + Omega
    S_1 <- MASS::ginv(S)

    # Kalman gain
    K <- Pstar %*% t(A) %*% S_1

    # Update P
    P  <- (Id_n - Pstar %*% t(A) %*%
             MASS::ginv(A %*% Pstar %*% t(A)) %*% A) %*% Pstar

    # Track change (unused, but kept for completeness)
    P.change <- P - old_P
    old_P <- P
  }

  # ----------------------------------------------------------------------------
  # 4. Build transition matrices for the full system (w_t, w_{t|t})
  # ----------------------------------------------------------------------------

  Theta_1 <- R %*% K %*% A %*% Phi_1
  Theta_2 <- R %*% ((Id_n - K %*% A) %*% Phi_1 + Lambda1_tilde)

  # Shocks mapping
  H11 <- (Id_n + Lambda0_tilde %*% R %*% K %*% A) %*% Sigma12
  H12 <- Lambda0_tilde %*% R %*% K %*% Omega12
  H21 <- R %*% K %*% A %*% Sigma12
  H22 <- R %*% K %*% Omega12

  # Save components back to model object
  model_sol <- model
  model_sol$R <- R
  model_sol$P <- P
  model_sol$K <- K
  model_sol$S <- S

  # ----------------------------------------------------------------------------
  # 5. Construct the joint VAR representation of (w_t, w_{t|t})
  # ----------------------------------------------------------------------------

  mu_ww <- rbind(mu_tilde, mu_tilde)

  Phi_ww <- rbind(cbind(Phi_1, Phi_2),
                  cbind(Theta_1, Theta_2))

  Sigma12_ww <- rbind(cbind(H11, H12),
                      cbind(H21, H22))

  Sigma_ww <- Sigma12_ww %*% t(Sigma12_ww)

  # ----------------------------------------------------------------------------
  # 6. Check stability (eigenvalues of the 2n × 2n transition matrix)
  # ----------------------------------------------------------------------------

  res_eig <- eigen(Phi_ww)
  if (sum(abs(res_eig$values) >= 1) > 0){
    warning("The resulting dynamics of w_t is not stationary
            (eigenvalues of Phi_tilde exceeding 1 in modulus).")
  }

  # ----------------------------------------------------------------------------
  # 7. Save final results
  # ----------------------------------------------------------------------------

  model_sol$mu_ww      <- mu_ww
  model_sol$Phi_ww     <- Phi_ww
  model_sol$Sigma_ww   <- Sigma_ww
  model_sol$Sigma12_ww <- Sigma12_ww

  return(model_sol)
}

