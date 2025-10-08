
compute_F_Shadow_Gaussian <-
  function(W,psi.parameterization,ell_bar,b,a,c,H){

  mu    <- psi.parameterization$mu
  Phi   <- psi.parameterization$Phi
  Sigma <- psi.parameterization$Sigma # Covariance matrix

  g <- function(x){
    return(x * pnorm(x) + dnorm(x))
  }

  n_x <- dim(Phi)[1] # number of states variables
  T <- dim(W)[1] # number of dates

  vecSS <- matrix(Sigma,ncol=1)

  # Initialization -------------------------------------------------------------
  all_b_n         <- matrix(0,H,1)
  all_b_bar_n     <- matrix(0,H,1)
  all_a_n         <- matrix(0,H,n_x)
  all_c_dot_n     <- matrix(0,H,1)
  all_c_n         <- matrix(0,H,n_x)
  all_sigma_n     <- matrix(0,H,1)

  all_a_tilde_n   <- matrix(0,H,n_x)
  all_c_tilde_n   <- matrix(0,H,n_x)

  # n = 1:
  Phi_low_n_1  <- matrix(0,n_x,n_x)
  Phi_exp_n_1  <- diag(n_x)
  sigma2_n_1 <- 0

  # Loop on maturities ---------------------------------------------------------

  for(i in 1:H){

    Phi_low_n <- Phi_low_n_1 + Phi_exp_n_1
    Phi_exp_n <- Phi %*% Phi_exp_n_1

    sigma2_n <- sigma2_n_1 + t(a) %*% Phi_exp_n_1 %*% Sigma %*% t(Phi_exp_n_1) %*% a

    b_bar_n <- b + t(a) %*% Phi_low_n %*% mu
    b_n     <- b_bar_n - .5 * t(a) %*% Phi_low_n %*% Sigma %*% t(Phi_low_n) %*% a
    c_dot_n <- t(c) %*% Phi_low_n %*% mu + .5 * t(c) %*% Phi_low_n %*% Sigma %*% t(Phi_low_n) %*% c

    # Store results:
    all_b_n[i]     <- b_n
    all_b_bar_n[i] <- b_bar_n
    all_a_n[i,]    <- t(Phi_exp_n) %*% a
    all_c_dot_n[i] <- c_dot_n
    all_c_n[i,]    <- t(Phi_exp_n) %*% c
    all_sigma_n[i] <- sqrt(sigma2_n)

    all_a_tilde_n[i,] <- t(Phi_low_n) %*% a
    all_c_tilde_n[i,] <- t(Phi_low_n) %*% c

    # For next iteration:
    Phi_low_n_1 <- Phi_low_n
    Phi_exp_n_1 <- Phi_exp_n
    sigma2_n_1  <- sigma2_n
  }

  vec1N <- matrix(1,H,1)
  vec1T <- matrix(1,T,1)
  vec1x <- matrix(1,n_x,1)

  aux_an_bnW    <- vec1T %*% t(all_b_n)     + W %*% t(all_a_n) - ell_bar
  aux_abarn_bnW <- vec1T %*% t(all_b_bar_n) + W %*% t(all_a_n) - ell_bar
  aux_s     <- aux_an_bnW    / (vec1T %*% t(all_sigma_n))
  aux_s_bar <- aux_abarn_bnW / (vec1T %*% t(all_sigma_n))

  aux_sigma <- vec1T %*% t(all_sigma_n)
  aux_pi    <- vec1T %*% t(all_c_dot_n) + W %*% t(all_c_n)

  aux_SS <- (all_a_tilde_n %x% t(vec1x)) * (t(vec1x) %x% all_c_tilde_n) * (vec1N %*% t(vecSS))
  aux_SS <- apply(aux_SS,1,sum)

  F_c_equal_0                <- ell_bar + aux_sigma * g(aux_s)
  indic_sigma0               <- which(all_sigma_n==0)
  F_c_equal_0[,indic_sigma0] <- ell_bar + pmax(aux_an_bnW[,indic_sigma0], 0)

  Phi_aux_s_bar <- pnorm(aux_s_bar)
  Phi_aux_s_bar[,indic_sigma0] <- 1 * (Phi_aux_s_bar[,indic_sigma0]>0)

  F_c <- F_c_equal_0 - aux_pi + Phi_aux_s_bar * (vec1T %*% t(aux_SS))

  Phi_aux_s <- pnorm(aux_s)
  Phi_aux_s[,indic_sigma0] <- 1 * (Phi_aux_s[,indic_sigma0]>0)

  phi_aux_s_bar <- dnorm(aux_s_bar) * (vec1T %*% t(aux_SS/all_sigma_n))
  phi_aux_s_bar[,indic_sigma0] <- 0

  all_nu_n_c_equal_0 <- Phi_aux_s
  all_nu_n_c         <- Phi_aux_s + phi_aux_s_bar

  return(list(F_c_equal_0 = F_c_equal_0,
              F = F_c
              # all_a_n     = all_a_n,
              # all_c_n     = all_c_n,
              # all_nu_n_c         = all_nu_n_c,
              # all_nu_n_c_equal_0 = all_nu_n_c_equal_0,
              # ptn = Phi_aux_s_bar,
              # all_b_bar_n = all_b_bar_n,
              # all_b_n = all_b_n,
              # all_sigma_n = all_sigma_n,
  ))
}

psi.QPoisson <- function(u,psi.parameterization){
  # Laplace transform of a process defined as follows:
  # N_t^+ ~ Poisson(x_t'·Omega^+·x_t) and N_t^- ~ Poisson(x_t'·Omega^-·x_t)
  # where Omega^+ and Omega^- are positive matrices and where
  # x_t follows a Gaussian VAR.
  Phi <- psi.parameterization$Phi
  Omega_plus  <- psi.parameterization$Omega_plus
  Omega_minus <- psi.parameterization$Omega_minus

  k <- dim(u)[2]
  n <- dim(Phi)[1]
  Mx <- make_Mnx(n)
  Kx <- make_Knx(n)

  u_aux <- matrix(NaN,n + n * (n + 1)/2,k)

  for (i in 1:k) {
    u_aux[1:n, i] <- matrix(u[1:n, i], ncol = 1) # this is v

    Vred <- u[(n + 1):(n + n * (n + 1)/2), i]
    vecV <- Mx %*% Vred
    V <- matrix(vecV, n, n)

    u_plus  <- u[n + n * (n + 1)/2 + 1, i]
    u_minus <- u[n + n * (n + 1)/2 + 2, i]

    V_aux <- V +
      (exp(u_plus) - 1)*Omega_plus +
      (exp(u_minus) - 1)*Omega_minus

    u_aux[(n + 1):(n + n * (n + 1)/2), i] <- Kx %*% matrix(V_aux,ncol=1)
  }

  res.QVAR <- psi.GaussianQVAR(u_aux, psi.parameterization)

  b <- res.QVAR$b
  a <- rbind(res.QVAR$a,
             u_plus,u_minus)

  return(list(a = a, b = b))
}

psi.VARG_Poisson <- function(u,psi.parameterization){
  # Laplace transform of a process defined as follows:
  # N_t^+ ~ Poisson(z_t^+) and N_t^- ~ Poisson(z_t^-),
  # where z_t = (z_t^+,z_t^-)' follows a bivariate VARG process.

  u <- matrix(u,nrow=4)
  k <- dim(u)[2]

  v <- matrix(u[1:2,],nrow=2)
  u_plus  <- u[3,]
  u_minus <- u[4,]

  # u_plus  <- 0
  # u_minus <- 0

  v_star <- v
  v_star[1,] <- v[1,] + exp(u_plus)  - 1
  v_star[2,] <- v[2,] + exp(u_minus) - 1

  res.VARG <- psi.VARG(v_star, psi.parameterization)

  b <- res.VARG$b
  a <- rbind(res.VARG$a,
             .99999*u_plus,.99999*u_minus)

  return(list(a = a, b = b))
}

psi.AR_CH <- function(u,psi.parameterization){
  # Laplace transform of a process defined as follows:
  # s_t = mu + phi·s_{t+1} + sqrt(z_t)·eps_t, with eps_t ~ N(0,1) and
  # z_t ARG(nu,beta,mu_z)

  mu  <- psi.parameterization$mu
  phi <- psi.parameterization$phi

  nu   <- psi.parameterization$nu
  mu_z <- psi.parameterization$mu_z
  beta <- psi.parameterization$beta

  u <- matrix(u,nrow=2)
  k <- dim(u)[2]

  u_s <- u[1,]
  u_z <- u[2,]

  u_z_tilde <- u_z + u_s^2/2

  b <- u_s*mu - nu*log(1 - mu_z*u_z_tilde)
  a <- matrix(NaN,2,k)
  a[1,] <- u_s * phi
  a[2,] <- beta*mu_z*u_z_tilde/(1 - mu_z*u_z_tilde)

  return(list(a = a, b = b))
}

varphi4G_SR_Gaussian <- function(x,parameterization,H){

  model   <- parameterization$model

  u <- parameterization$u
  v <- parameterization$v # determine thresholds (dimension nw x 1)
  i.v.x <- matrix(v,ncol=1) %*% matrix(c(1i*x),nrow=1)

  nw <- dim(i.v.x)[1]
  q  <- dim(i.v.x)[2]

  u2 <- matrix(0, nw, q)
  u1 <- matrix(u, nw, q) + i.v.x

  res_reverse <- reverse.MHLT(psi.GaussianVAR,
                              u1 = u1,
                              u2 = u2,
                              H = H,
                              psi.parameterization = model)
  A <- res_reverse$A
  B <- res_reverse$B
  return(list(A=A,B=B))
}

varphi4G_SR_QPoisson <- function(x,parameterization,H){

  model   <- parameterization$model

  u <- parameterization$u
  v <- parameterization$v # determine thresholds (dimension nw x 1)
  i.v.x <- matrix(v,ncol=1) %*% matrix(c(1i*x),nrow=1)

  nw <- dim(i.v.x)[1]
  q  <- dim(i.v.x)[2]

  u2 <- matrix(0, nw, q)
  u1 <- matrix(u, nw, q) + i.v.x

  res_reverse <- reverse.MHLT(psi.PoissonSR,
                              u1 = u1,
                              u2 = u2,
                              H = H,
                              psi.parameterization = model)
  A <- res_reverse$A
  B <- res_reverse$B
  return(list(A=A,B=B))
}

compute_F_Shadow_affine <- function(W,psi,psi.parameterization,
                                    ell_bar=0,
                                    b,a,c,chi=-1,
                                    H,
                                    eps = 10^(-6), # to compute dG
                                    max_x = 2000,
                                    dx_statio = 10,
                                    min_dx = 1e-06,
                                    nb_x1 = 1000,
                                    du = 1e-06 # To numerically compute Gamma0
                                    #            and Gamma1 (condit. variance)
){
  # W is of dimension TT x n, where TT is the number of dates.
  # psi.parameterization contains the model specifications

  TT <- dim(W)[1] # number of dates.
  vec1TT <- matrix(1,TT,1)

  n_w <- dim(W)[2]
  n_a <- length(a)
  n_c <- length(c)
  if((n_a != n_w)|(n_c != n_w)){
    print("Error: a and c should be vectors of dimenson n_w x 1, and W a #dates x n_w matrix.")
    return(1)
  }

  varphi <- function(x,parameterization,H){
    # This function is used in truncated.payof

    u <- parameterization$u
    v <- parameterization$v # determine thresholds (dimension nw x 1)
    i.v.x <- matrix(v,ncol=1) %*% matrix(c(1i*x),nrow=1)

    nw <- dim(i.v.x)[1]
    q  <- dim(i.v.x)[2]

    u2 <- matrix(0, nw, q)
    u1 <- matrix(u, nw, q) + i.v.x

    res_reverse <- reverse.MHLT(psi,
                                u1 = u1,
                                u2 = u2,
                                H = H,
                                psi.parameterization =
                                  parameterization$model)
    A <- res_reverse$A
    B <- res_reverse$B
    return(list(A=A,B=B))
  }

  # Compute G0 and dG: ---------------------------------------------------------
  parameterization <- list(
    model = psi.parameterization,
    u = matrix(0*a,ncol=1),
    v = matrix(1*a,ncol=1)
  )
  res_truncated0 <- truncated.payoff(W,b.matrix = matrix(ell_bar-b,H,1),
                                     varphi = varphi,
                                     parameterization = parameterization,
                                     max_x = max_x,
                                     dx_statio = dx_statio,
                                     min_dx = min_dx,
                                     nb_x1 = nb_x1)
  parameterization$u <- matrix(eps*a,ncol=1)
  res_truncatedeps <- truncated.payoff(W,b.matrix = matrix(ell_bar-b,H,1),
                                       varphi = varphi,
                                       parameterization = parameterization,
                                       max_x = max_x,
                                       dx_statio = dx_statio,
                                       min_dx = min_dx,
                                       nb_x1 = nb_x1)
  G0   <- matrix(res_truncated0,nrow=TT)
  Geps <- matrix(res_truncatedeps,nrow=TT)
  dG <- (Geps - G0)/eps
  # ----------------------------------------------------------------------------

  psi.parameterization$n_w <- n_w # this is required in 'compute_expect_variance'
  res_EV <- compute_expect_variance(psi,psi.parameterization,du = du)

  # Compute E_t(a'w_{t+n}): ----------------------------------------------------
  res_EV_0a <- compute_expect_variance_H(res_EV,
                                         theta1 = t(a),theta2 = t(0*a),H = H)
  C0 <- matrix(res_EV_0a$all_C0,nrow=1)
  C1 <- matrix(res_EV_0a$all_C1,ncol=H)
  E_aW <- matrix(1,TT,1) %*% C0 + W %*% C1

  # Compute E_t(c'w_{t+n}): ----------------------------------------------------
  res_EV_0c <- compute_expect_variance_H(res_EV,
                                         theta1 = t(c),theta2 = t(0*c),H = H)
  C0 <- matrix(res_EV_0c$all_C0,nrow=1)
  C1 <- matrix(res_EV_0c$all_C1,ncol=H)
  E_cW <- matrix(1,TT,1) %*% C0 + W %*% C1

  # Compute sigma2_n(c-a,w_{t}): -----------------------------------------------
  res_EV_c_a <- compute_expect_variance_H(res_EV,
                                          theta1 = t(c+chi*a),theta2 = t(c+chi*a),H = H)
  sigma2_c_a <- vec1TT %*% matrix(res_EV_c_a$all_Gamma0,nrow=1) +
    W %*% matrix(res_EV_c_a$all_Gamma1,n_w,H)

  # Compute sigma2_n(c,w_{t}): -----------------------------------------------
  res_EV_cc <- compute_expect_variance_H(res_EV,
                                         theta1 = t(c),theta2 = t(c),H = H)
  sigma2_cc <- vec1TT %*% matrix(res_EV_cc$all_Gamma0,nrow=1) +
    W %*% matrix(res_EV_cc$all_Gamma1,n_w,H)

  # Compute sigma2_n(a,w_{t}): -----------------------------------------------
  res_EV_aa <- compute_expect_variance_H(res_EV,
                                         theta1 = t(a),theta2 = t(a),H = H)
  sigma2_aa <- vec1TT %*% matrix(res_EV_aa$all_Gamma0,nrow=1) +
    W %*% matrix(res_EV_aa$all_Gamma1,n_w,H)

  delta_sigma2_c_a <- sigma2_c_a - cbind(0,sigma2_c_a[,1:(H-1)])
  delta_sigma2_cc  <- sigma2_cc  - cbind(0,sigma2_cc[,1:(H-1)])
  delta_sigma2_aa  <- sigma2_aa  - cbind(0,sigma2_aa[,1:(H-1)])

  # Combine all sub-results to compute F:
  F_c <- b*(1 - G0) + E_aW - E_cW - dG +
    (-1/2*(1-G0)*delta_sigma2_c_a) +
    (-1/2*G0*delta_sigma2_cc)

  F_c_equal_0 <- - chi*b*(1 - G0) - chi*E_aW + chi*dG +
    (-1/2*(1-G0)*delta_sigma2_aa)

  return(list(F_c_equal_0 = F_c_equal_0,
              F = F_c))
}


