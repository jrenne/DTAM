
#---------------------------------------------------------------------------
# Simulate data from a Switching State-Space (SSRS) model
#   Reference: Kim and Nelson (1999), Section 5.1
#
# Model structure:
#   Measurement equation:
#       y_t = H_{s_t} β_t + A_{s_t} Z_t + e_t,   e_t ~ N(0, R_{s_t})
#
#   State transition equation:
#       β_t = μ_{s_t} + F_{s_t} β_{t-1} + G_{s_t} v_t,   v_t ~ N(0, Q_{s_t})
#
# where the regime s_t ∈ {1, …, M} follows a first-order Markov chain
# with transition matrix P (rows summing to one).
#---------------------------------------------------------------------------

simul_SSRS <- function(SSmodel, s0 = NaN, beta0, Z) {

  #---------------------------------------------------------------------------
  # 1. Extract model dimensions and components
  #---------------------------------------------------------------------------

  T <- dim(Z)[1]         # Number of time periods (rows in Z)
  H  <- SSmodel$H        # Measurement matrices (N × J × M)
  A  <- SSmodel$A        # Coefficients on exogenous regressors (N × K × M)
  R  <- SSmodel$R        # Covariance of measurement shocks (N × N × M)
  mu <- SSmodel$mu       # Regime-dependent intercepts in the state eq (J × 1 × M)
  F  <- SSmodel$F        # State transition matrices (J × J × M)
  G  <- SSmodel$G        # Loadings on state shocks (J × L × M)
  Q  <- SSmodel$Q        # Covariance of state shocks (L × L × M)
  P  <- SSmodel$P        # Regime transition matrix (M × M, rows sum to 1)

  N <- dim(H)[1]         # Number of observed variables
  J <- dim(F)[1]         # Number of latent state variables
  K <- dim(A)[2]         # Number of exogenous variables
  L <- dim(G)[2]         # Number of shocks driving β_t
  M <- dim(H)[3]         # Number of regimes

  #---------------------------------------------------------------------------
  # 2. Simulate the regime sequence {s_t}
  #    Each s_t is a 1×M one-hot vector (e.g., [1,0,0]) indicating regime.
  #---------------------------------------------------------------------------

  s.simul <- simul_RS(P, TT = T, ini_state = s0)

  #---------------------------------------------------------------------------
  # 3. Initialize storage matrices for simulated series
  #---------------------------------------------------------------------------

  y.simul    <- matrix(NaN, T, N)  # Simulated observations
  beta.simul <- matrix(NaN, T, J)  # Simulated latent states

  beta <- beta0  # Initial state vector β₀

  #---------------------------------------------------------------------------
  # 4. Simulation loop over time
  #---------------------------------------------------------------------------

  for (t in 1:T) {

    # Identify the current regime index s ∈ {1,…,M}
    s <- which(s.simul[t, ] == 1)

    #----------------------------------------
    # 4.1 State transition equation
    #       β_t = μ_s + F_s β_{t-1} + G_s v_t
    #----------------------------------------

    v <- t(chol(Q[,,s])) %*% rnorm(L)     # Draw v_t ~ N(0, Q_s)
    beta <- mu[, 1, s] + F[,,s] %*% beta + G[,,s] %*% v

    #----------------------------------------
    # 4.2 Measurement equation
    #       y_t = H_s β_t + A_s Z_t + e_t
    #----------------------------------------

    e <- t(chol(R[,,s])) %*% rnorm(N)     # Draw e_t ~ N(0, R_s)
    y <- matrix(H[,,s],N,J) %*% beta + A[,,s] %*% matrix(Z[t,], ncol = 1) + e

    #----------------------------------------
    # 4.3 Store simulated values
    #----------------------------------------

    y.simul[t, ]    <- y
    beta.simul[t, ] <- beta
  }

  #---------------------------------------------------------------------------
  # 5. Return results
  #---------------------------------------------------------------------------

  return(list(
    s.simul    = s.simul,    # Simulated regime indicators (T × M)
    y.simul    = y.simul,    # Simulated observables (T × N)
    beta.simul = beta.simul  # Simulated latent states (T × J)
  ))
}


#---------------------------------------------------------------------------
# Kim94(): Kim's (1994) Collapsing Filter for Markov-Switching State-Space Models
#
# References:
#   - Kim, C.-J. (1994). "Dynamic Linear Models with Markov-Switching."
#     Journal of Econometrics, 60(1–2), 1–22.
#   - Kim and Nelson (1999), Section 5.1.
#
# Model:
#   Measurement: y_t = H_{s_t} β_t + A_{s_t} Z_t + e_t,   e_t ~ N(0, R_{s_t})
#   Transition:  β_t = μ_{s_t} + F_{s_t} β_{t−1} + G_{s_t} v_t,   v_t ~ N(0, Q_{s_t})
#
# where s_t ∈ {1, …, M} follows a first-order Markov chain with transition
# probabilities matrix P (each row sums to 1).
#
# This function performs filtering using:
#   - The Hamilton filter for regime probabilities.
#   - A Kalman filter for continuous states within each regime pair (i,j).
#   - Kim’s (1994) "collapse" step to keep mixture size tractable.
#
# Output:
#   - loglik: log-likelihood of the sample
#   - all.beta.tt: filtered state means E[β_t | Y_t]
#   - all.P.tt: filtered state variances vec(V[β_t | Y_t])
#   - all.xi.t: filtered regime probabilities P(s_t | Y_t)
#---------------------------------------------------------------------------

Kim94 <- function(SSmodel, Y.matrix, Z.matrix, beta0, P0, xi0 = NaN) {

  #---------------------------------------------------------------------------
  # 1. Helper: Multivariate normal density function (for regime pair likelihoods)
  #---------------------------------------------------------------------------
  N.multi.pdf <- function(eta, Omega) {
    # eta: N × k matrix of residuals
    # Omega: N × N × k array of covariance matrices
    k <- dim(eta)[2]
    N <- dim(eta)[1]
    res <- apply(matrix(1:k, k, 1), 1, function(j) {
      if(dim(Omega)[1]>1){
        detOmegaj <- det(Omega[,, j])
      }else{
        detOmegaj <- abs(Omega[,, j])
      }
      -0.5 * t(eta[, j]) %*% solve(Omega[,, j]) %*% eta[, j]-
        .5*log((2 * pi)^N * detOmegaj)
    })
    return(res)
  }

  #---------------------------------------------------------------------------
  # 2. Extract model parameters and dimensions
  #---------------------------------------------------------------------------

  H  <- SSmodel$H   # Measurement matrices (N × J × M)
  A  <- SSmodel$A   # Exogenous regressor matrices (N × K × M)
  R  <- SSmodel$R   # Measurement error covariances (N × N × M)
  mu <- SSmodel$mu  # State intercepts (J × 1 × M)
  F  <- SSmodel$F   # State transition matrices (J × J × M)
  G  <- SSmodel$G   # Shock loadings (J × L × M)
  Q  <- SSmodel$Q   # State shock covariances (L × L × M)
  P  <- SSmodel$P   # Regime transition matrix (M × M, rows sum to 1)

  N <- dim(H)[1]    # Number of observed variables
  J <- dim(F)[1]    # Number of state variables
  K <- dim(A)[2]    # Number of exogenous regressors
  L <- dim(G)[2]    # Number of shocks v_t
  M <- dim(H)[3]    # Number of regimes
  T <- dim(Y.matrix)[1]  # Sample size

  regimes <- matrix(1:M, ncol = 1)  # Convenience matrix of regime indices

  #---------------------------------------------------------------------------
  # 3. Initialize regime-dependent containers
  #---------------------------------------------------------------------------

  beta_i.tt     <- matrix(NaN, J, M)        # Mean of β_t | Y_t, for each regime i
  P_i.tt        <- matrix(NaN, J^2, M)      # Covariance (vectorized)
  f.ij          <- matrix(NaN, M, M)        # Densities f(y_t | s_{t-1}=i, s_t=j)
  beta_ij.tt_1  <- array(NaN, c(J, M, M))   # Predicted state means (i → j)
  beta_ij.tt    <- array(NaN, c(J, M, M))   # Updated state means (i → j)
  eta_ij.tt_1   <- array(NaN, c(N, M, M))   # Forecast errors (i → j)
  P_ij.tt_1     <- array(NaN, c(J^2, M, M)) # Predicted covariances (i → j)
  P_ij.tt       <- array(NaN, c(J^2, M, M)) # Updated covariances (i → j)
  f_ij.tt_1     <- array(NaN, c(N^2, M, M)) # Prediction error covariances (i → j)

  #---------------------------------------------------------------------------
  # 4. Storage for outputs
  #---------------------------------------------------------------------------

  all.beta.tt <- matrix(NaN, T, J)   # Filtered state means (aggregated)
  all.P.tt    <- matrix(NaN, T, J^2) # Filtered covariances (aggregated)
  all.xi.t    <- matrix(NaN, T, M)   # Filtered regime probabilities

  loglik <- 0                        # Initialize log-likelihood accumulator

  #---------------------------------------------------------------------------
  # 5. Initialization of state and covariance
  #---------------------------------------------------------------------------

  beta_i.tt <- matrix(beta0, J, M)   # Initial state mean for all regimes
  P_i.tt    <- matrix(P0, J^2, M)    # Initial covariance for all regimes

  #---------------------------------------------------------------------------
  # 6. Initialization of regime probabilities
  #---------------------------------------------------------------------------

  if (!is.na(xi0[1])) {
    xi.t <- xi0
  } else {
    # Compute stationary distribution by repeated squaring of P'
    P_2h <- t(P)
    for (i in 1:10) {
      P_2h <- P_2h %*% P_2h
    }
    xi.t <- P_2h[, 1]
  }

  #---------------------------------------------------------------------------
  # 7. Main filtering recursion
  #---------------------------------------------------------------------------

  for (t in 1:T) {

    #--- Save previous step results ------------------------------------------
    beta_i.t_1t_1 <- beta_i.tt
    P_i.t_1t_1    <- P_i.tt
    xi.t_1        <- xi.t

    y_t <- matrix(Y.matrix[t, ], ncol = 1)
    z_t <- matrix(Z.matrix[t, ], ncol = 1)

    # Joint prior probabilities P(s_t=j, s_{t-1}=i | Y_{t−1})
    Pr.St.St_1.psi_t_1 <- (xi.t_1 %*% matrix(1, 1, M)) * P

    #--- Regime-by-regime forecasting and updating ---------------------------
    for (i in 1:M) {

      # Eqs. (5.10)–(5.13): prediction
      beta_ij.tt_1[,, i] <- apply(regimes, 1, function(j) {
        mu[,, j] + F[,, j] %*% beta_i.t_1t_1[, i]
      })

      P_ij.tt_1[,, i] <- apply(regimes, 1, function(j) {
        F[,, j] %*% matrix(P_i.t_1t_1[, i], J, J) %*% t(F[,, j]) +
          G[,, j] %*% Q[,, j] %*% t(G[,, j])
      })

      eta_ij.tt_1[,, i] <- apply(regimes, 1, function(j) {
        y_t - matrix(H[,, j],N,J) %*% beta_ij.tt_1[, j, i] - A[,, j] %*% z_t
      })

      f_ij.tt_1[,, i] <- apply(regimes, 1, function(j) {
        matrix(H[,, j],N,J) %*% matrix(P_ij.tt_1[, j, i], J, J) %*%
          t(matrix(H[,, j],N,J)) + R[,, j]
      })

      # Conditional densities f(y_t | s_{t−1}=i, s_t=j)
      f.ij[i, ] <- N.multi.pdf(matrix(eta_ij.tt_1[,, i],N, M), R)

      # Eqs. (5.14)–(5.15): updating
      beta_ij.tt[,, i] <- beta_ij.tt_1[,, i] +
        apply(regimes, 1, function(j) {
          matrix(P_ij.tt_1[, j, i], J, J) %*% t(matrix(H[,, j],N,J)) %*%
            solve(matrix(f_ij.tt_1[, j, i], N, N)) %*% eta_ij.tt_1[, j, i]
        })

      P_ij.tt[,, i] <- apply(regimes, 1, function(j) {
        (diag(J) -
           matrix(P_ij.tt_1[, j, i], J, J) %*%
           t(matrix(H[,, j],N,J)) %*% solve(matrix(f_ij.tt_1[, j, i], N, N)) %*%
           matrix(H[,, j],N,J)) %*% matrix(P_ij.tt_1[, j, i], J, J)
      })
    }

    #--- Hamilton filter update for regime probabilities ---------------------
    f.y_t.psi_t_1 <- sum(exp(f.ij) * Pr.St.St_1.psi_t_1)
    Pr.St.St_1.psi_t <- exp(f.ij) * Pr.St.St_1.psi_t_1 / f.y_t.psi_t_1
    xi.t <- apply(Pr.St.St_1.psi_t, 2, sum)

    #--- Kim (1994) collapsing step ------------------------------------------
    # Collapse mixture components by matching first two conditional moments

    beta_i.tt <- apply(regimes, 1, function(j) {
      beta_ij.tt[, j, ] %*% Pr.St.St_1.psi_t[, j]
    })

    beta_i.tt.Array <- array(beta_i.tt, c(J, M, M))
    diff.beta_ij.tt <- beta_i.tt.Array - beta_ij.tt
    beta.beta.ij.tt <- apply(diff.beta_ij.tt, c(2, 3), function(x) { x %x% x })

    P_i.tt <- apply(regimes, 1, function(j) {
      (P_ij.tt[, j, ] + beta.beta.ij.tt[, j, ]) %*% Pr.St.St_1.psi_t[, j]
    })

    #--- Store results --------------------------------------------------------
    all.beta.tt[t, ] <- beta_i.tt %*% xi.t
    all.P.tt[t, ]    <- P_i.tt %*% xi.t
    all.xi.t[t, ]    <- xi.t

    #--- Update log-likelihood -----------------------------------------------
    loglik <- loglik + log(f.y_t.psi_t_1)
  }

  #---------------------------------------------------------------------------
  # 8. Return filtered results
  #---------------------------------------------------------------------------

  return(list(
    loglik      = loglik,
    all.beta.tt = all.beta.tt,
    all.P.tt    = all.P.tt,
    all.xi.t    = all.xi.t
  ))
}


# N <- 2 # number of observed variables
# J <- 2 # number of latent variables
# L <- J # number of shocks v
# K <- 1 # number of exogenous variables
# M <- 3 # number of regimes
#
# SSmodel <- list()
# SSmodel$H <- array(rnorm(N*J*M),c(N,J,M))
# SSmodel$A <- array(rnorm(N*K*M),c(N,K,M))
# SSmodel$R  <- array(2*diag(N),c(N,N,M))
#
# SSmodel$mu <- array(rnorm(J*M),c(J,1,M))
# SSmodel$F  <- array(.8*diag(J),c(J,J,M))
# SSmodel$G  <- array(diag(J),c(J,J,M))
# SSmodel$Q  <- array(.1*diag(J),c(J,J,M))
# SSmodel$P  <- t(matrix(c(.8,.15,.05,
#                          .05,.9,.05,
#                          0,.1,.9),M,M)) # matrix of transition probabilities, row components sum to 1.
#
#
# beta0 <- matrix(0,J,1)
# T <- 200
# Z <- matrix(1,T,K)
#
# res.simul <- simul_SSRS(SSmodel=SSmodel,
#                         beta0=beta0,Z=Z)
#
# # eta <- matrix(rnorm(8),4,2)
# # Omega <- array(diag(4),c(4,4,2))
# # N.multi.pdf(eta,Omega)
#
# res.filter <- Kim94(SSmodel,Y.matrix=res.simul$y.simul,Z,beta0,P0=diag(J))
#
# par(mfrow=c(2,2))
# plot(res.simul$s.simul[,1],type="l")
# lines(res.filter$all.xi.t[,1],col="red")
# plot(res.simul$s.simul[,2],type="l")
# lines(res.filter$all.xi.t[,2],col="red")
# plot(res.simul$beta.simul[,1],type="l")
# lines(res.filter$all.beta.tt[,1],col="red")
# plot(res.simul$beta.simul[,2],type="l")
# lines(res.filter$all.beta.tt[,2],col="red")



