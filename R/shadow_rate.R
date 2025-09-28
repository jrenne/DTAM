# simul.var <- function(Model,nb.sim,x0=NaN){
#   n <- dim(Model$rho)[1]
#   if(is.na(x0[1])){
#     x0 <- solve(diag(n) - Model$rho) %*% Model$mu
#   }
#
#   X <- c(x0)
#   x <- x0
#   for(t in 2:nb.sim){
#     x <- Model$mu + Model$rho %*% x + Model$Sigma %*% rnorm(n)
#     X <- rbind(X,c(x))
#   }
#   return(X)
# }
#
# compute.price.WX <- function(Model,X,max.H){
#   # max.H is the maximum considered maturity, expressed in number of periods.
#
#   SS <- Model$Sigma %*% t(Model$Sigma)
#
#   n <- dim(Model$rho)[1]
#
#   # For h=1:
#   sum.rho.j_1 <- diag(n)
#   rho.j <- diag(n)
#
#   sigma.h.2 <- 0
#
#   all.a.bar <- NULL
#   all.b <- NULL
#   all.a <- NULL
#   all.sigma <- NULL
#
#   for(h in 1:max.H){
#     sigma.h.2 <- sigma.h.2 + matrix(Model$delta.1,nrow=1) %*% rho.j %*% SS %*%
#       t(rho.j) %*% matrix(Model$delta.1,ncol=1)
#
#     a.h.bar <- Model$delta.0 + matrix(Model$delta.1,nrow=1) %*% sum.rho.j_1 %*% Model$mu
#     a.h <- a.h.bar - .5*matrix(Model$delta.1,nrow=1) %*% sum.rho.j_1 %*%
#       SS %*% t(sum.rho.j_1) %*% matrix(Model$delta.1,ncol=1)
#
#     rho.j <- rho.j %*% Model$rho
#
#     b.h <- matrix(Model$delta.1,nrow=1) %*% rho.j
#
#     all.a.bar <- c(all.a.bar,a.h.bar)
#     all.b <- rbind(all.b, b.h)
#     all.a <- c(all.a, a.h)
#     all.sigma <- c(all.sigma,sqrt(sigma.h.2))
#
#     sum.rho.j_1 <- sum.rho.j_1 + rho.j # for next iteration
#   }
#
#   # Computation of forward rates:
#
#   vec.1 <- matrix(1,dim(X)[1],1)
#
#   aux <- vec.1 %*% matrix(all.a,nrow=1) +  X %*% t(all.b) - Model$r.bar
#   aux <- aux / (vec.1 %*% matrix(all.sigma,nrow=1))
#
#   ggg <- function(z){
#     return(z * pnorm(z) + dnorm(z))
#   }
#
#   vec.f <- Model$r.bar + (vec.1 %*% matrix(all.sigma,nrow=1)) * ggg(aux)
#
#   vec.y <- t(apply(vec.f,1,cumsum))
#   vec.y <- vec.y / (vec.1 %*% matrix(1:max.H,nrow=1))
#
#   return(list(
#     all.a.bar = all.a.bar,
#     all.a = all.a,
#     all.b = all.b,
#     all.sigma = all.sigma,
#     aux = aux,
#     vec.f = vec.f,
#     vec.y = vec.y
#   ))
# }


apply_Lemma <- function(mu,Phi,Sigma,i_bar,b,a,c,X,max_h){

  g <- function(x){
    return(x * pnorm(x) + dnorm(x))
  }

  n_x <- dim(Phi)[1] # number of states variables
  T <- dim(X)[1] # number of dates

  SS    <- Sigma %*% t(Sigma)
  vecSS <- matrix(SS,ncol=1)

  # Initialization -------------------------------------------------------------
  all_b_n         <- matrix(0,max_h,1)
  all_b_bar_n     <- matrix(0,max_h,1)
  all_a_n         <- matrix(0,max_h,n_x)
  all_c_dot_n     <- matrix(0,max_h,1)
  all_c_n         <- matrix(0,max_h,n_x)
  all_sigma_n     <- matrix(0,max_h,1)

  all_a_tilde_n   <- matrix(0,max_h,n_x)
  all_c_tilde_n   <- matrix(0,max_h,n_x)

  # n = 1:
  Phi_low_n_1  <- matrix(0,n_x,n_x)
  Phi_exp_n_1  <- diag(n_x)
  sigma2_n_1 <- 0

  # Loop on maturities ---------------------------------------------------------

  for(i in 1:max_h){

    Phi_low_n <- Phi_low_n_1 + Phi_exp_n_1
    Phi_exp_n <- Phi %*% Phi_exp_n_1

    sigma2_n <- sigma2_n_1 + t(a) %*% Phi_exp_n_1 %*% SS %*% t(Phi_exp_n_1) %*% a

    b_bar_n <- b + t(a) %*% Phi_low_n %*% mu
    b_n     <- b_bar_n - .5 * t(a) %*% Phi_low_n %*% SS %*% t(Phi_low_n) %*% a
    c_dot_n <- t(c) %*% Phi_low_n %*% mu + .5 * t(c) %*% Phi_low_n %*% SS %*% t(Phi_low_n) %*% c

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

  vec1N <- matrix(1,max_h,1)
  vec1T <- matrix(1,T,1)
  vec1x <- matrix(1,n_x,1)

  aux_an_bnX    <- vec1T %*% t(all_b_n)     + X %*% t(all_a_n) - i_bar
  aux_abarn_bnX <- vec1T %*% t(all_b_bar_n) + X %*% t(all_a_n) - i_bar
  aux_s     <- aux_an_bnX    / (vec1T %*% t(all_sigma_n))
  aux_s_bar <- aux_abarn_bnX / (vec1T %*% t(all_sigma_n))

  aux_sigma <- vec1T %*% t(all_sigma_n)
  aux_pi    <- vec1T %*% t(all_c_dot_n) + X %*% t(all_c_n)

  aux_SS <- (all_a_tilde_n %x% t(vec1x)) * (t(vec1x) %x% all_c_tilde_n) * (vec1N %*% t(vecSS))
  aux_SS <- apply(aux_SS,1,sum)

  F_c_equal_0                <- i_bar + aux_sigma * g(aux_s)
  indic_sigma0               <- which(all_sigma_n==0)
  F_c_equal_0[,indic_sigma0] <- i_bar + pmax(aux_an_bnX[,indic_sigma0], 0)

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
              F_c         = F_c,
              all_a_n     = all_a_n,
              all_c_n     = all_c_n,
              all_nu_n_c         = all_nu_n_c,
              all_nu_n_c_equal_0 = all_nu_n_c_equal_0,
              ptn = Phi_aux_s_bar,
              all_b_bar_n = all_b_bar_n,
              all_b_n = all_b_n,
              all_sigma_n = all_sigma_n,
              all_a_n = all_a_n))
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

  k <- dim(u)[2]

  v <- matrix(u[1:2,],nrow=2)
  u_plus  <- u[3,]
  u_minus <- u[4,]

  v_star <- v
  v_star[1,] <- v[1,] + exp(u_plus)  - 1
  v_star[2,] <- v[2,] + exp(u_minus) - 1

  res.VARG <- psi.VARG(v_star, psi.parameterization)

  b <- res.VARG$b
  a <- rbind(res.VARG$a,
             u_plus,u_minus)

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

varphi4G_SR <- function(x,psi,parameterization,H){
  # This employs the (single-horizon) Laplace transform and
  #    makes it suitable to be used by function compute_F_Shadow_affine.
  # 'parameterization' is a list containing:
  #     - a 'model' (second argument of the 'psi' function)
  #     - XXXX

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
                                parameterization$psi.parameterization)
  A <- res_reverse$A
  B <- res_reverse$B
  return(list(A=A,B=B))
}

compute_F_Shadow_affine <- function(W,psi,psi.parameterization,
                                    ell_bar,b,a,c,
                                    H,
                                    eps = 10^(-6), # to compute dG
                                    max_x = 2000,
                                    dx_statio = 1,
                                    min_dx = 1e-06,
                                    nb_x1 = 1000){
  # W is of dimension TT x n, where TT is the number of dates.

  TT <- dim(W)[1] # number of dates.

  # Compute G0 and dG:

  parameterization <- list(
    psi.parameterization = psi.parameterization,
    u = matrix(0*xi1,ncol=1),
    v = matrix(1*xi1,ncol=1)
  )
  res_truncated0 <- truncated.payoff(W,b.matrix = matrix(i_bar-xi0,H,1),
                                     varphi = varphi4G_SR_Gaussian,
                                     parameterization = parameterization,
                                     max_x = max_x,
                                     dx_statio = dx_statio,
                                     min_dx = min_dx,
                                     nb_x1 = nb_x1)
  parameterization$u <- matrix(eps*xi1,ncol=1)
  res_truncatedeps <- truncated.payoff(W,b.matrix = matrix(i_bar-xi0,H,1),
                                       varphi = varphi4G_SR_Gaussian,
                                       parameterization = parameterization,
                                       max_x = max_x,
                                       dx_statio = dx_statio,
                                       min_dx = min_dx,
                                       nb_x1 = nb_x1)
  G0   <- matrix(res_truncated0,nrow=TT)
  Geps <- matrix(res_truncatedeps,nrow=TT)

  dG <- (Geps - G0)/eps

  model$n_w <- 2
  model <- XXXX
  res_EV <- compute_expect_variance_H(compute_expect_variance(psi,model),
                                      theta1 = t(xi1),theta2 = t(0*xi1),H = H)

  C0 <- matrix(res_EV$all_C0,nrow=1)
  C1 <- matrix(res_EV$all_C1,ncol=H)
  E_aW <- matrix(1,TT,1) %*% C0 + W %*% C1

  res_EV2 <- compute_expect_variance_H(compute_expect_variance(psi,model),
                                       theta1 = t(xi1),theta2 = t(xi1),H = H)
  delta_sigma2 <- c(res_EV2$all_Gamma0) - c(0,c(res_EV2$all_Gamma0)[1:(H-1)])

  F <- xi0*(1 - G0) + E_aW - dG - 1/2*(1-G0)*t(matrix(delta_sigma2,H,TT))

  return(F)
}


