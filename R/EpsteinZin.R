

solve_EZ_SDF <- function(model,psi,Ew=NaN,z_bar_ini=5,
                         nb_loop_z_bar=20,nb_loop_mu_z1=100){
  # This procedure computes the SDF in the CES-CRRA Epstein-Zin case.
  # -- Inputs ------------------------------------------------------------------
  # model is a list containing the model parameters, including:
  # - the parameterization of preferences, i.e., rho, gamma, delta.
  # - the parameterization of consumption growth, i.e., mu_c0 and mu_c1
  # - the dimension of w, parameter: n_w (used to compute Ew numerically)
  # psi is the log Laplace transform of the state vector w_t; it can itself use
  #     model as an argument.
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
  if(is.na(Ew[1])){
    du   <- 10^(-8)
    mu <- matrix(0,n_w,1)
    Phi <- matrix(0,n_w,n_w)
    for(i in 1:n_w){
      u_du <- matrix(0,model$n_w,1)
      u_du[i] <- du
      psi_w <- psi(u_du,model)
      mu[i] <- psi_w$b/du
      Phi[i,] <- psi_w$a/du
    }
    Ew <- solve(diag(n_w) - Phi) %*% mu
  }

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

  return(
    list(Ew = Ew, Phi = Phi, mu = mu,
         mu_z0 = mu_z0, mu_z1 = mu_z1, z_bar = z_bar,
         eta0 = eta0, eta1 = eta1,
         alpha = alpha,
         kappa0 = kappa0, kappa1 = kappa1))
}

psi.BansalShaliastovich <- function (u,model){
  # This function computes the log Laplace transform associated with
  # the model by Bansal and Shaliastovich (2012, RFS)
  n <- dim(model$Pi)[1]
  q <- dim(model$nu)[1]
  K <- dim(u)[2]
  u.x <- matrix(u[1:n, ], nrow = n)
  u.z <- matrix(u[(n + 1):(n + q), ], nrow = q)
  u.epsilon <- matrix(u[(n + q + 1):(n + q + 3), ], nrow = 3)
  b <- 0.5 * t(u.x * u.x) %*% matrix(model$chi.0, ncol = 1) +
    0.5 * t(u.z * u.z) %*% matrix(model$sigma.w^2, ncol = 1) +
    0.5 * ((t(u.epsilon) %*% model$Omega) * t(u.epsilon)) %*%
    matrix(1, 3, 1)
  a <- rbind(t(model$Pi) %*% u.x, 0.5 * t(model$nu) %*%
               model$chi.1 %*% (u.x * u.x) + t(model$nu) %*% u.z, matrix(0,
                                                                         3, K))
  return(list(a = a,b = b))
}

