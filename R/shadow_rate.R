simul.var <- function(Model,nb.sim,x0=NaN){
  n <- dim(Model$rho)[1]
  if(is.na(x0[1])){
    x0 <- solve(diag(n) - Model$rho) %*% Model$mu
  }

  X <- c(x0)
  x <- x0
  for(t in 2:nb.sim){
    x <- Model$mu + Model$rho %*% x + Model$Sigma %*% rnorm(n)
    X <- rbind(X,c(x))
  }
  return(X)
}


compute.price.WX <- function(Model,X,max.H){
  # max.H is the maximum considered maturity, expressed in number of periods.

  SS <- Model$Sigma %*% t(Model$Sigma)

  n <- dim(Model$rho)[1]

  # For h=1:
  sum.rho.j_1 <- diag(n)
  rho.j <- diag(n)

  sigma.h.2 <- 0

  all.a.bar <- NULL
  all.b <- NULL
  all.a <- NULL
  all.sigma <- NULL

  for(h in 1:max.H){
    sigma.h.2 <- sigma.h.2 + matrix(Model$delta.1,nrow=1) %*% rho.j %*% SS %*%
      t(rho.j) %*% matrix(Model$delta.1,ncol=1)

    a.h.bar <- Model$delta.0 + matrix(Model$delta.1,nrow=1) %*% sum.rho.j_1 %*% Model$mu
    a.h <- a.h.bar - .5*matrix(Model$delta.1,nrow=1) %*% sum.rho.j_1 %*%
      SS %*% t(sum.rho.j_1) %*% matrix(Model$delta.1,ncol=1)

    rho.j <- rho.j %*% Model$rho

    b.h <- matrix(Model$delta.1,nrow=1) %*% rho.j

    all.a.bar <- c(all.a.bar,a.h.bar)
    all.b <- rbind(all.b, b.h)
    all.a <- c(all.a, a.h)
    all.sigma <- c(all.sigma,sqrt(sigma.h.2))

    sum.rho.j_1 <- sum.rho.j_1 + rho.j # for next iteration
  }

  # Computation of forward rates:

  vec.1 <- matrix(1,dim(X)[1],1)

  aux <- vec.1 %*% matrix(all.a,nrow=1) +  X %*% t(all.b) - Model$r.bar
  aux <- aux / (vec.1 %*% matrix(all.sigma,nrow=1))

  ggg <- function(z){
    return(z * pnorm(z) + dnorm(z))
  }

  vec.f <- Model$r.bar + (vec.1 %*% matrix(all.sigma,nrow=1)) * ggg(aux)

  vec.y <- t(apply(vec.f,1,cumsum))
  vec.y <- vec.y / (vec.1 %*% matrix(1:max.H,nrow=1))

  return(list(
    all.a.bar = all.a.bar,
    all.a = all.a,
    all.b = all.b,
    all.sigma = all.sigma,
    aux = aux,
    vec.f = vec.f,
    vec.y = vec.y
  ))
}


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
              all_a_n = all_a_n))
}


# mu <- matrix(c(.01,0,0),ncol=1)
# Phi <- .8 * diag(3)
# Phi[2,1] <- .1
# Sigma12 <- .02*diag(3)
# Sigma12[3,1] <- .0
# Sigma <- Sigma12 %*% t(Sigma12)
#
# TT <- 200
# Model <- list(mu=mu,Phi=Phi,Sigma = Sigma)
# X <- simul.GVAR(Model,nb.sim=TT)
#
# b <- .0
# a <- matrix(c(1,0,0),ncol=1)
# c <- matrix(c(0,0,1),ncol=1)
#
# EX <- X %*% t(Phi) + matrix(1,TT,1) %*% t(mu)
#
# res <- apply_Lemma(mu,Phi,Sigma,i_bar=0,b,a,c,X,max_h=10)
#
# r <- b + X %*% a
# plot(r,type="l")
# lines(res$F_c_equal_0[,1],col="red")
# lines(b + EX %*% a,col="blue",lty=2)


psi.PoissonSR <- function(u,psi.parameterization){
  mu <- matrix(psi.parameterization$mu, ncol = 1)
  Phi <- psi.parameterization$Phi
  Sigma <- psi.parameterization$Sigma # Covariance matrix
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


# rho <- .98
# mu <- .16
#
# mu <- matrix( mu*(1-rho),2,1)
# Phi <- rho*diag(2)
# Sigma12 <- .02*matrix(c(1,-1),2,1)
# Sigma <- Sigma12 %*% t(Sigma12)
#
# model <- list(mu=mu,Phi=Phi,Sigma=Sigma,Sigma12=Sigma12,n_w=2)
#
# TT = 200
#
# simX <- simul.GVAR(model,nb.sim = TT)
#
# par(mfrow=c(2,2))
# par(plt=c(.1,.95,.1,.9))
# plot(simX[,1],type="l",ylim=c(min(0,min(simX)),max(simX)))
# lines(simX[,2],type="l",col="red")
# abline(h=0)
#
# lambda_plus  <- simX[,1]^2
# lambda_minus <- simX[,2]^2
#
# plot(lambda_plus,type="l",ylim=c(0,max(simX^2)))
# lines(lambda_minus,type="l",col="red")
# abline(h=0)
#
# DeltaNplus <- rpois(n = TT,lambda_plus)
# DeltaNminus <- rpois(n = TT,lambda_plus)
#
# Nplus  <- cumsum(DeltaNplus)
# Nminus <- cumsum(DeltaNminus)
#
# W <- cbind(simX,
#            ((simX %x% matrix(1,1,2))*(matrix(1,1,2) %x% simX))[,c(1,2,4)],
#            Nplus,Nminus)
#
# plot(Nplus,type="l")
# lines(Nminus,col="red")
#
# r <- 1 + .25*(Nplus - Nminus)
# plot(r,type="l",lwd=2)
# abline(h=0)
#
#
# n <- dim(model$Phi)[1]
#
# xi0 <- .01
# xi1 <- matrix(c(rep(0,n),
#                 rep(0,n*(n+1)/2),
#                 .25/100,-.25/100),ncol=1)
# H <- 10
#
# model$Omega_plus  <- diag(c(1,0))
# model$Omega_minus <- diag(c(0,1))
#
# u <- - xi1
# psi.PoissonSR(u,model)
#
# res <- compute_AB_classical(xi0, xi1, kappa0 = NaN, kappa1 = NaN,
#                             H=H, psi=psi.PoissonSR,
#                             psi.parameterization=model)
#
# yields <- matrix(1,TT,1) %*% matrix(res$b,nrow=1) +
#   W %*% matrix(res$a,dim(res$a)[1],dim(res$a)[3])
#
# t <- 40
# lines((t+1):(t+H),100*yields[t,],col="red",lwd=2)
#
#
# parameterization <- list(
#   model = model,
#   H = H,
#   u = matrix(0*xi1,ncol=1),
#   v = matrix(1*xi1,ncol=1)
# )
#
# #varphi4G_SR(1,parameterization)
#
# # res_truncated <- truncated.payoff(W,b.matrix = matrix(-xi0,H,1),H = H,
# #                                   varphi = varphi4G_SR,
# #                                   parameterization = parameterization)
#
# i_bar <- 0
#
# rho <- .98
# mu <- .02
# mu <- matrix( mu*(1-rho),2,1)
# Phi <- rho*diag(2)
# Sigma12 <- .015*matrix(c(1,-1),2,1)
# Sigma <- Sigma12 %*% t(Sigma12)
#
# model <- list(mu=mu,Phi=Phi,Sigma=Sigma,Sigma12=Sigma12,n_w=2)
#
# TT = 200
# simX <- simul.GVAR(model,nb.sim = TT)
#
# par(mfrow=c(2,1))
# abline(h=0,col="grey")
# plot(simX[,1],type="l")
#
# xi1 <- matrix(c(1,0),ncol=1)
# xi0 <- 0
#
# parameterization <- list(
#   model = model,
#   H = H,
#   u = matrix(0*xi1,ncol=1),
#   v = matrix(1*xi1,ncol=1)
# )
# res_truncated0 <- truncated.payoff(simX,b.matrix = matrix(i_bar-xi0,H,1),
#                                    varphi = varphi4G_SR_Gaussian,
#                                    parameterization = parameterization,
#                                    max_x = 2000,
#                                    dx_statio = 1,
#                                    min_dx = 1e-06,
#                                    nb_x1 = 1000)
# eps <- 10^(-6)
# parameterization$u <- matrix(eps*xi1,ncol=1)
# res_truncatedeps <- truncated.payoff(simX,b.matrix = matrix(i_bar-xi0,H,1),
#                                      varphi = varphi4G_SR_Gaussian,
#                                      parameterization = parameterization,
#                                      max_x = 2000,
#                                      dx_statio = 1,
#                                      min_dx = 1e-06,
#                                      nb_x1 = 1000)
# G0   <- matrix(res_truncated0,nrow=TT)
# Geps <- matrix(res_truncatedeps,nrow=TT)
#
# dG <- (Geps - G0)/eps
#
# h <- 10
# plot(G0[,h],type="l",ylim=c(0,1))
# lines(G0[,1],col="red")
#
# model$n_w <- 2
# res_EV <- compute_expect_variance_H(compute_expect_variance(psi.GaussianVAR,model),
#                                     theta1 = t(xi1),theta2 = t(0*xi1),H = 10)
#
# C0 <- matrix(res_EV$all_C0,nrow=1)
# C1 <- matrix(res_EV$all_C1,ncol=H)
# E_aW <- matrix(1,TT,1) %*% C0 + simX %*% C1
#
# res_EV2 <- compute_expect_variance_H(compute_expect_variance(psi.GaussianVAR,model),
#                                      theta1 = t(xi1),theta2 = t(xi1),H = 10)
# delta_sigma2 <- c(res_EV2$all_Gamma0) - c(0,c(res_EV2$all_Gamma0)[1:(H-1)])
#
# F <- xi0*(1 - G0) + E_aW - dG - 1/2*(1-G0)*t(matrix(delta_sigma2,H,TT))
#
# h <- 10
#
# plot(F[,h],type="l")
#
# # Compare with Extended Wu-Xia:
# b <- xi0
# a <- xi1
# c <- 0*xi1
# res_GaussianLemma <- apply_Lemma(model$mu,
#                                  model$Phi,
#                                  model$Sigma12,i_bar=i_bar,b,a,c,simX,max_h=H)
# lines(res_GaussianLemma$F_c_equal_0[,h],col="red")
#
#
# h <- 10
# plot(simX[,1],type="l")
# plot(1-res_GaussianLemma$ptn[,h],type="l")
# lines(G0[,h],col="red")
#
# sigma <- t(sqrt(matrix(res_EV$all_Gamma0,H,TT)))
# mu <- E_aW
# beta <- - mu/sigma
# EaW1_aW <- pnorm(beta)*mu - dnorm(beta)*sigma
#
# h <- 10
# plot(EaW1_aW[,h],type="l")
# lines(dG[,h],col="red")
#
#
#
