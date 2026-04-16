# #---------------------------------------------------------------------------
# # Build SSmodel object for trend–cycle model with drift in trend
# # and regime-dependent cycle (compatible with Kim94 and simul_SSRS)
# #---------------------------------------------------------------------------
#
# build_SSmodel_trend_cycle_drift <- function(mu_tau, phi, sigma_tau, sigma_c,
#                                             mu_c1, mu_c2, pi11, pi22) {
#
#   #---------------------------------------------------------------------------
#   # 1. Dimensions
#   #---------------------------------------------------------------------------
#   N <- 2  # number of observed variables (y_t)
#   J <- 2  # number of state variables (tau_t, c_t)
#   K <- 1  # number of exogenous variables (variable will be = 0)
#   L <- J  # number of shocks (ε_tau,t, ε_c,t)
#   M <- 2  # number of regimes
#
#   #---------------------------------------------------------------------------
#   # 2. Measurement equation: y_t = [1 1] * [tau_t, c_t]'
#   #---------------------------------------------------------------------------
#   H <- array(0, c(N, J, M))
#   for (m in 1:M) H[1,,m] <- matrix(c(1, 1), nrow = 1)
#   for (m in 1:M) H[2,,m] <- matrix(c(1, 0), nrow = 1)
#
#   # No exogenous regressors
#   A <- array(1, c(N, K, M))
#
#   # Measurement noise (set small if negligible)
#   R <- array(0, c(N, N, M))
#   for (m in 1:M) R[,,m] <- 1e-3*diag(N)
#
#   #---------------------------------------------------------------------------
#   # 3. State transition equations
#   # beta_t = mu_s + F_s * beta_{t-1} + G_s * v_t
#   #---------------------------------------------------------------------------
#   mu <- array(0, c(J, 1, M))
#   mu[, , 1] <- matrix(c(mu_tau, mu_c1), ncol = 1)   # regime 1 (recession)
#   mu[, , 2] <- matrix(c(mu_tau, mu_c2), ncol = 1)   # regime 2 (expansion)
#
#   # Transition matrix F (same across regimes)
#   F <- array(0, c(J, J, M))
#   for (m in 1:M) {
#     F[,,m] <- matrix(c(1, 0,
#                        0, phi), nrow = 2, byrow = TRUE)
#   }
#
#   # Shock loading matrix (identity)
#   G <- array(0, c(J, L, M))
#   for (m in 1:M) G[,,m] <- diag(c(sigma_tau, sigma_c))
#
#   # Covariance of state shocks
#   Q <- array(0, c(L, L, M))
#   for (m in 1:M) {
#     Q[,,m] <- diag(L)
#   }
#
#   #---------------------------------------------------------------------------
#   # 4. Regime transition matrix
#   #---------------------------------------------------------------------------
#   P <- matrix(c(pi11, 1 - pi11,
#                 1 - pi22, pi22), nrow = 2, byrow = TRUE)
#
#   #---------------------------------------------------------------------------
#   # 5. Return SSmodel list
#   #---------------------------------------------------------------------------
#   SSmodel <- list(H = H, A = A, R = R, mu = mu, F = F, G = G, Q = Q, P = P)
#   return(SSmodel)
# }
#
# # Parameters
# mu_tau    <- 0.005      # drift of trend
# phi       <- 0.7      # cycle persistence
# sigma_tau <- 0.01
# sigma_c   <- 0.1
# pi11      <- 0.8
# pi22      <- 0.95
#
# mu_c1     <- -0.01 # regime 1 (recession)
# mu_c2     <- (-mu_c1)*(pi11/pi22) # regime 2 (expansion)
#
# # mu_c1 <- 0
# # mu_c2 <- mu_c1
#
#
# # Build model
# SSmodel <- build_SSmodel_trend_cycle_drift(mu_tau, phi, sigma_tau, sigma_c,
#                                            mu_c1, mu_c2, pi11, pi22)
#
# # Simulate and filter
# T <- 200
# Z <- matrix(0, T, 1)
# beta0 <- matrix(0, 2, 1)
#
# sim <- simul_SSRS(SSmodel, beta0 = beta0, Z = Z)
#
# res <- Kim94(SSmodel,
#              Y.matrix = sim$y.simul,
#              Z.matrix = Z,
#              beta0=sim$y.simul[1],
#              P0 = diag(2))
#
# par(mfrow=c(2,2))
# par(plt=c(.1,.95,.1,.9))
# plot(sim$y.simul[,1],type="l")
# lines(res$all.beta.tt[,1],col="red")
# plot(sim$beta.simul[,2],type="l")
# lines(res$all.beta.tt[,2],col="red")
#
# plot(sim$s.simul[,1],type="l")
# lines(res$all.xi.t[,1],col="red")
#
