

solve_EZ_SDF <- function(model,psi,Ew=NaN,z_bar_ini=5,
                         nb_loop_z_bar=20,nb_loop_mu_z1=100){
  # This procedure computes the SDF in the CES-CRRA Epstein-Zin case.
  # -- Inputs ------------------------------------------------------------------
  # model is a list containing the model parameters, including:
  # - the parameterization of preferences, i.e., rho, gamma, delta.
  # - the parameterization of consumption growth, i.e., mu_c0 and mu_c1
  # - the dimension of w, parameter: n_w (used to compute Ew numerically)
  # psi is the log Laplace transform of the state vector w_t; it can itself use
  #     model as an argument (i.e., model contains the process parameterization).
  # Ew is the mean of w_t. If NaN, it is computed numerically.
  # ----------------------------------------------------------------------------

  rho   <- model$rho
  gamma <- model$gamma
  delta <- model$delta

  mu_c0 <- model$mu_c0
  mu_c1 <- model$mu_c1

  theta <- (1 - gamma)/(1 - rho)

  n_w  <- model$n_w

  # Computation of Ew, using the Laplace Transform (Ew = dPsi(u)/du at u=0) ----
  res_moments <- compute_expect_variance(psi,model)
  Ew  <- res_moments$Ew
  Vw  <- res_moments$Vw
  Phi <- res_moments$Phi
  mu  <- res_moments$mu

  # Loop on z_bar: -------------------------------------------------------------
  z_bar <- z_bar_ini

  for(iii in 1:nb_loop_z_bar){
    kappa1 <- c(exp(z_bar)/(1 + exp(z_bar)))
    kappa0 <- c(log(1 + exp(z_bar)) - kappa1 * z_bar)

    # Solving for mu_z1 (fixed-point problem recursively solved) ---------------
    # Initial value based on VAR representation of w_t process:
    mu_z1 <- (1 - rho) * solve(diag(n_w) - kappa1*t(Phi)) %*% t(Phi) %*% mu_c1
    for(jjj in 1:nb_loop_mu_z1){
      u <- matrix(theta*kappa1*mu_z1 - theta*(rho-1)*mu_c1,ncol=1)
      psi_w <- psi(u,model)
      mu_z1 <- 1/theta * psi_w$a
      mu_z0 <- 1/(1-kappa1)*(log(delta) - (rho-1)*mu_c0 + kappa0 +
                               1/theta*psi_w$b)
    }
    # update z_bar:
    z_bar <- mu_z0 + t(mu_z1) %*% Ew
  }

  # Determine vector of prices of risk: ----------------------------------------
  alpha <- (theta-1)*kappa1*mu_z1 - (theta*rho-theta+1)*mu_c1

  # Determine real short-term rate specification: ------------------------------
  psi_w_alpha <- psi(alpha,model)
  eta0 <- - theta*log(delta) + (theta*rho-theta+1)*mu_c0 -
    (theta-1)*(kappa0 + kappa1*mu_z0 - mu_z0) - psi_w_alpha$b
  eta1 <- (theta-1)*mu_z1 - psi_w_alpha$a

  model_solved <- model

  model_solved$Ew     <- Ew
  model_solved$Vw     <- Vw
  model_solved$Phi    <- Phi
  model_solved$mu     <- mu
  model_solved$mu_z0  <- c(mu_z0)
  model_solved$mu_z1  <- mu_z1
  model_solved$z_bar  <- z_bar
  model_solved$eta0   <- c(eta0)
  model_solved$eta1   <- eta1
  model_solved$alpha  <- alpha
  model_solved$kappa0 <- kappa0
  model_solved$kappa1 <- kappa1

  return(model_solved)
}

# psi.BansalShaliastovich <- function (u,model){
#   # This function computes the log Laplace transform associated with
#   # the model by Bansal and Shaliastovich (2012, RFS)
#   n <- dim(model$Pi)[1]
#   q <- dim(model$nu)[1]
#   K <- dim(u)[2]
#   u.x <- matrix(u[1:n, ], nrow = n)
#   u.z <- matrix(u[(n + 1):(n + q), ], nrow = q)
#   u.epsilon <- matrix(u[(n + q + 1):(n + q + 3), ], nrow = 3)
#   b <- 0.5 * t(u.x * u.x) %*% matrix(model$chi.0, ncol = 1) +
#     0.5 * t(u.z * u.z) %*% matrix(model$sigma.w^2, ncol = 1) +
#     0.5 * ((t(u.epsilon) %*% model$Omega) * t(u.epsilon)) %*%
#     matrix(1, 3, 1)
#   a <- rbind(t(model$Pi) %*% u.x, 0.5 * t(model$nu) %*%
#                model$chi.1 %*% (u.x * u.x) + t(model$nu) %*% u.z, matrix(0,
#                                                                          3, K))
#   return(list(a = a,b = b))
# }



solve_EZ_stock_return <- function(model,psi,Ew=NaN,z_s_bar_ini=5,
                                  nb_loop_z_bar=20,nb_loop_mu_z1=100){
  # This procedure computes an affine approximation to stock return in
  # the context of Epstein-Zin preferences.
  # The approach is that of Bansal and Yaron (2004), based on the linearization
  # of returns à la Campbell and Shiller (1988).
  # -- Inputs ------------------------------------------------------------------
  # model_solved is a list containing the model parameters, including:
  # - the parameterization of preferences, i.e., rho, gamma, delta
  # - the parameterization of consumption growth, i.e., mu_c0 and mu_c1
  # - the parameterization of dividend growth, i.e., mu_d0 and mu_d1.
  # - the dimension of w, parameter: n_w (used to compute Ew numerically)
  # In addition, model_solved contains the SDF specification, that includes:
  # - the vector of prices of risk (alpha)
  # - the specification of the short-term rate (scalar eta0 and vector eta1).
  # psi is the log Laplace transform of the state vector w_t; it can itself use
  #     model as an argument (i.e., model contains the process parameterization).
  # Ew is the mean of w_t. If NaN, it is computed numerically.
  # ----------------------------------------------------------------------------

  n_w <- model$n_w
  Phi <- model$Phi
  Ew  <- model$Ew
  Phi <- model$Phi

  mu_d0 <- model$mu_d0
  mu_d1 <- model$mu_d1

  eta0  <- model$eta0
  eta1  <- model$eta1
  alpha <- model$alpha

  # Loop on z_bar: -------------------------------------------------------------
  z_s_bar <- z_s_bar_ini

  for(iii in 1:nb_loop_z_bar){
    kappa1s <- c(exp(z_s_bar)/(1 + exp(z_s_bar)))
    kappa0s <- c(log(1 + exp(z_s_bar)) - kappa1s * z_s_bar)

    # Solving for mu_z1 (fixed-point problem recursively solved) ---------------
    # Initial value based on VAR representation of w_t process:
    mu_s1 = solve(diag(n_w) - kappa1s*t(Phi)) %*% (t(Phi) %*% mu_d1 - eta1)
    for(jjj in 1:nb_loop_mu_z1){
      u <- matrix(alpha + kappa1s*mu_s1 + mu_d1,ncol=1)
      psi_w_alpha_plus <- psi(u,model)
      psi_w_alpha      <- psi(alpha,model)
      mu_s1 <- psi_w_alpha_plus$a - eta1 - psi_w_alpha$a
      mu_s0 <- 1/(1-kappa1s)*(psi_w_alpha_plus$b - eta0 - psi_w_alpha$b +
                                kappa0s + mu_d0)
    }
    # update z_bar:
    z_s_bar <- mu_s0 + t(mu_s1) %*% Ew
  }

  # Deduce specification of r_{s,t+1}:
  mu_rs0 <- kappa0s + kappa1s*mu_s0 - mu_s0 + mu_d0
  mu_rs1 <- - mu_s1
  mu_rs2 <- kappa1s * mu_s1 + mu_d1

  model_solved <- model

  model_solved$mu_rs0 <- c(mu_rs0)
  model_solved$mu_rs1 <- mu_rs1
  model_solved$mu_rs2 <- mu_rs2

  model_solved$mu_s0 <- c(mu_s0)
  model_solved$mu_s1 <- mu_s1

  model_solved$kappa0s <- kappa0s
  model_solved$kappa1s <- kappa1s

  model_solved$z_s_bar <- z_s_bar

  return(model_solved)
}


psi.BansalShaliastovich <- function (u,model){
  # This function computes the log Laplace transform associated with
  # the model by Bansal and Shaliastovich (2012, RFS)
  n <- dim(model$Pi)[1]
  q <- dim(model$nu)[1]
  u.x <- matrix(u[1:n, ], nrow = n)
  u.z <- matrix(u[(n + 1):(n + q), ], nrow = q)
  b <- 0.5 * t(u.x * u.x) %*% matrix(model$chi.0, ncol = 1) +
    0.5 * t(u.z * u.z) %*% matrix(model$sigma.w^2, ncol = 1)
  a <- rbind(t(model$Pi) %*% u.x, 0.5 *
               model$chi.1 %*% (u.x * u.x) + t(model$nu) %*% u.z)
  return(list(a = a,b = b))
}

