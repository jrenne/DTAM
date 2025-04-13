
reverse.MHLT <- function(psi,u1,u2=NaN,H,psi.parameterization){
  # compute the multi-horizon Laplace transform in the reverse-order case
  # That is, we consider:
  # E_t(exp(u_2'w_{t+1}+...+u_2'w_{t+h-1}+u_1'w_{t+h}))
  # Inputs: psi is the (one-period) conditional Laplace transform of process w_t
  # The function can evaluate the MHLT in parallel for k different vectors
  #    u1 and u2 (i.e., u1 and u2, as inputs can be of dimension n x k with k > 1)
  if(is.na(u2[1])){# in that case, u2 <- u1
    u2 <- u1
  }
  A = NULL
  B = NULL
  A.h_1 <- 0
  B.h_1 <- 0
  n <- dim(u1)[1]
  k <- dim(u1)[2]
  A <- array(NaN,c(n,k,H))
  B <- array(NaN,c(1,k,H))
  for(i in 1:H){
    if(i==1){u <- u1}else{u <- u2}
    psi.u <- psi(u + A.h_1,psi.parameterization)
    A.h <- psi.u$a
    B.h <- psi.u$b + B.h_1
    # save results
    A[,,i] <- A.h
    B[,,i] <- B.h
    # for next iteration
    A.h_1 <- A.h
    B.h_1 <- B.h
  }
  matrix1overH <- -1/array((1:H) %x% rep(1,n*k),c(n,k,H))
  a <- A*matrix1overH
  matrix1overH <- -1/array((1:H) %x% rep(1,1*k),c(1,k,H))
  b <- B*matrix1overH
  return(list(A=A,B=B,a=a,b=b))
}

psi.GaussianVAR <- function(u,psi.parameterization){
  # Laplace transform of a Gaussian VAR:
  # w_t = mu + Phi w_{t-1} + epsilon_{t}, where epsilon_{t}~N(0,Sigma)
  # If w_t is n-dimensional (of dimension n, say), u is of dimension n x k
  #    (i.e., we can compute k LT in parallel)
  # WARNING: Sigma is a covariance matrix.
  mu    <- psi.parameterization$mu
  Phi   <- psi.parameterization$Phi
  Sigma <- psi.parameterization$Sigma
  a <- t(Phi) %*% u
  b <- t(u) %*% mu + .5 * t(apply(u,2,function(x){x %x% x})) %*% c(Sigma)
  return(list(a=a,b=b))
}
# # Check:
# mu <- matrix(1:3,ncol=1)
# Phi <- .9 * diag(3)
# Sigma <- diag(3)
# Sigma[1,3] <- .5
# Sigma[3,1] <- .5
# psi.parameterization <- list(
#   mu = mu,
#   Phi = Phi,
#   Sigma = Sigma
# )
# u <- matrix(1:12,nrow=3)
# psi.GaussianVAR(u,psi.parameterization)
# # Check reverse-order multi-horizon LT:
# u1 <- matrix(1,3,2)
# u2 <- u1/2
# reverse.MHLT(psi.GaussianVAR,u1,u2,H,psi.parameterization)




simul.ARG <- function(nb.sim,mu,nu,rho,alpha=0,w0=NaN){
  # This function simulates an ARG or an ARG0 process

  unc.mean <- (alpha + nu) * mu/(1-rho)

  if(is.na(w0)){
    w_0 <- unc.mean
  }else{
    w_0 <- w0
  }
  W <- w_0
  w <- w_0
  for(t in 2:nb.sim){
    z <- rpois(1,rho*w/mu + alpha)
    w <- mu * rgamma(1,shape=nu+z)
    W <- c(W,w)
  }
  return(W)
}

simul.compound.poisson <- function(nb.sim,Gamma,Pi,lambda,w0=NaN){
  # This function simulates a compound Poisson process

  if(is.na(w0)){
    w_0 <- 1 * Gamma
  }else{
    w_0 <- w0
  }

  W <- w_0
  w <- w_0
  for(t in 2:nb.sim){
    z <- sum(rbernoulli(n = w/Gamma, p = Pi))
    eps <- rpois(1,lambda)
    w <- Gamma * (z + eps)
    W <- c(W,w)
  }
  return(W)
}

rbernoulli <- function(n,p){
  return(1*(runif(n)<p))
}

simul.MS.AR <- function(nb.sim,mu.1,mu.2,rho.1,rho.2,sigma.1,sigma.2,P,w0=NaN){
  # This function simulates a Markov-Switching AR process
  # s is valued in {1,2}
  # p(1,2) = P( s_{t}=2 | s_{t-1}=1 )
  max.2.power <- 8
  Pk <- t(P)
  for(i in 1:max.2.power){
    Pk <- Pk %*% Pk
  }
  if(Pk[1,1]<.5){
    s <- 2
  }else{
    s <- 1
  }

  S <- s

  if(is.na(w0)){
    w_0 <- Pk[1,1] * mu.1/(1-rho.1) + Pk[2,1] * mu.2/(1-rho.2)
  }else{
    w_0 <- w0
  }

  W <- w_0
  w <- w_0

  for(t in 2:nb.sim){
    u <- runif(1)
    if(s==1){
      if(u>P[1,1]){
        s <- 2
      }
    }else{
      if(u<P[2,1]){
        s <- 1
      }
    }
    S <- c(S,s)

    if(s==1){
      w <- mu.1 + rho.1 * w + sigma.1 * rnorm(1)
    }else{
      w <- mu.2 + rho.2 * w + sigma.2 * rnorm(1)
    }

    W <- c(W,w)
  }
  return(list(S=S,W=W,Pk=Pk))
}

make.pdf <- function(model,values.of.variable,
                     psi,
                     max_x = 10000,
                     min_x = 10^(-8),nb_x=5000){
  # psi is a function calculating the Laplace transform (LT) of a univariate random variable
  # psi admits 2 arguments: model and u. model is a list of parameters; u is the
  # vector of points where one wants to evaluate the LT.

  x <- matrix(exp(seq(log(min_x),log(max_x),length.out=nb_x)),ncol=1)
  gamma <- matrix(values.of.variable,nrow=1)

  cdf.values <- Fourier.psi(model,gamma,x,psi)

  nb.values.variable <- length(values.of.variable)

  scale.variable.values <- seq(c(values.of.variable)[1],
                               tail(c(values.of.variable),1),
                               length.out = nb.values.variable+1)

  cdf <- pmax(pmin(cdf.values,1),0)
  for(i in 2:length(cdf)){
    if(cdf[i]<cdf[i-1]){
      cdf[i] <- cdf[i-1]
    }
  }

  tmp <- splinefun(x=values.of.variable, y=cdf, method="hyman")

  fitted.cdf.values <- tmp(scale.variable.values)
  fitted.pdf.values <- diff(fitted.cdf.values)*mean(diff(values.of.variable,1))

  return(fitted.pdf.values)
}

Fourier.psi <- function(model,gamma,x,
                        psi){
  # x has to be a column vector
  # gamma has to be a row vector
  u <- c(1i*x)
  psi_eval <- matrix(psi(model,u),length(x),length(gamma))
  dx <- matrix(x-c(0,x[1:length(x)-1]),length(x),length(gamma))
  fx <- Im(psi_eval*exp(-1i*x %*% gamma))/matrix(x,length(x),length(gamma))*dx

  f  <- 1/2 - 1/pi * apply(fx,2,sum)

  f.cdf <- pmax(pmin(f,1),0)
  for(i in 2:length(f)){
    if(f.cdf[i]<f.cdf[i-1]){
      f.cdf[i] <- f.cdf[i-1]}}

  return(f.cdf)
}

simul.VAR <- function (Model, nb.sim, x0 = NaN)
{
  n <- dim(Model$Phi)[1]
  if (is.na(x0[1])) {
    x0 <- solve(diag(n) - Model$Phi) %*% Model$mu
  }
  X <- c(x0)
  x <- x0
  for (t in 2:nb.sim) {
    x <- Model$mu + Model$Phi %*% x + Model$Sigma %*% rnorm(n)
    X <- rbind(X, c(x))
  }
  return(X)
}

compute_expect_variance <- function(psi,model){
  # This function computes the conditional and unconditional expectations and
  # varances of an affine process.
  # ----------------------------------------------------------------------------
  # The conditional mean is given by:
  # E_t(w_{t+1}) = mu + Phi.w_t
  # The conditional variance is given by:
  # vec[Var_t(w_{t+1})] = Gamma_0 + Gamma_1.w_t
  # ----------------------------------------------------------------------------

  n_w  <- model$n_w

  # Computation of Ew, using the Laplace Transform (Ew = dPsi(u)/du at u=0) ----
  du   <- 10^(-9)
  mu <- matrix(0,n_w,1)
  Phi <- matrix(0,n_w,n_w)
  Gamma0 <- matrix(0,n_w*n_w,1)
  Gamma1 <- matrix(0,n_w*n_w,n_w)

  for(i in 1:n_w){
    u_plus_i <- matrix(0,model$n_w,1)
    u_plus_i[i] <- du
    u_minu_i <- matrix(0,model$n_w,1)
    u_minu_i[i] <- -du
    psi_plus_i <- psi(u_plus_i,model)
    psi_minu_i <- psi(u_minu_i,model)
    # Compute specification of first-order conditional moment:
    mu[i]   <- .5 * (psi_plus_i$b - psi_minu_i$b)/du
    Phi[i,] <- .5 * (psi_plus_i$a - psi_minu_i$a)/du
    for(j in 1:n_w){
      u_plus_i_plus_j <- u_plus_i
      u_plus_i_minu_j <- u_plus_i
      u_plus_i_plus_j[j] <- u_plus_i_plus_j[j] + du
      u_plus_i_minu_j[j] <- u_plus_i_minu_j[j] - du
      u_minu_i_plus_j <- u_minu_i
      u_minu_i_minu_j <- u_minu_i
      u_minu_i_plus_j[j] <- u_minu_i_plus_j[j] + du
      u_minu_i_minu_j[j] <- u_minu_i_minu_j[j] - du
      psi_plus_i_plus_j <- psi(u_plus_i_plus_j,model)
      psi_plus_i_minu_j <- psi(u_plus_i_minu_j,model)
      psi_minu_i_plus_j <- psi(u_minu_i_plus_j,model)
      psi_minu_i_minu_j <- psi(u_minu_i_minu_j,model)
      # Compute specification of second order moment:
      Gamma0[n_w*(i-1)+j,1] <-
        (psi_plus_i_plus_j$b - psi_plus_i_minu_j$b - psi_minu_i_plus_j$b + psi_minu_i_minu_j$b)/
        (4*du^2)
      Gamma1[n_w*(i-1)+j,]  <-
        (psi_plus_i_plus_j$a - psi_plus_i_minu_j$a - psi_minu_i_plus_j$a + psi_minu_i_minu_j$a)/
        (4*du^2)
    }
    # Compute unconditional expectation:
    Ew <- solve(diag(n_w) - Phi) %*% mu
    # Compute unconditional variance:
    Vw <- matrix(solve(diag(n_w^2) - Phi %x% Phi) %*% (Gamma0 + Gamma1 %*% Ew),
                 n_w,n_w)
  }

  return(list(
    mu = mu,
    Phi = Phi,
    Gamma0 = Gamma0,
    Gamma1 = Gamma1,
    Ew = Ew,
    Vw = Vw
  ))
}

# # Check
# modelVAR <- list(mu = matrix(c(1,0),2,1),
#                  Phi = matrix(c(.6,.2,-.3,.8),2,2),
#                  Sigma = diag(2),
#                  n_w=2)
# compute_expect_variance(psi.GaussianVAR,modelVAR)
# solve(diag(4) - modelVAR$Phi %x% modelVAR$Phi) %*% c(modelVAR$Sigma)


psi.VARG <- function(u,model){
  # Laplace transform of a VARG with conditionally independent components.
  # Conditionally on w_t,w_{t-1},..., the distribution of the jth component is:
  #  w_{j,t} ~ gamma_{nu_j}(alpha_j + beta_j'w_{t-1},mu_j),
  # In model, the vectors alpha, nu, mu gather the n values
  #  of the associated parameters. The rows of matrix beta are the beta_j'
  # If w_t is n-dimensional (of dimension n, say), u is of dimension n x k
  #    (i.e., we can compute k LT in parallel)
  alpha <- model$alpha
  beta  <- model$beta
  nu    <- model$nu
  mu    <- model$mu
  n <- dim(u)[1]
  k <- dim(u)[2]
  mu_matrix <- matrix(mu,n,k)
  a <- t(beta) %*% ((u * mu_matrix)/(1 - u * mu_matrix))
  b <- t((u * mu_matrix)/(1 - u * mu_matrix)) %*% alpha -
    t(log(1 - u * mu_matrix)) %*% nu
  return(list(a=a,b=b))
}

compute_P_bar_VARG <- function(model,gamma,H=10,indic_delta1 = NaN){
  # gamma has to be a vector with the same dimension as w_t
  # model contains, in particular, the parameterization of the VARG model.
  # Prices are computed directly for all defaultable entities.
  # If indic_delta1 is an integer, it indicates the component of w_t that
  #  corresponds to delta_{e=1}. If NaN, delta_{e=1} is supposed to correspond
  #  to the first non-zero entry of vector nu.
  xi0 <- model$xi0
  xi1 <- model$xi1

  alpha <- model$alpha
  nu    <- model$nu
  mu    <- model$mu
  beta  <- model$beta

  if(is.na(indic_delta1[1])){
    indic_delta1 <- which(nu==0)[1]
  }
  n.w <- dim(beta)[1]
  n.e <- n.w - indic_delta1 + 1
  n.y <- n.w - n.e

  e_tilde <- rbind(matrix(0,n.y,n.e),
                   diag(n.e))
  gamma_matrix <- matrix(gamma,n.w,n.e)
  xi1_matrix   <- matrix(xi1,n.w,n.e)

  u <- -10000

  u1 <- gamma_matrix + u*e_tilde
  u2 <- - xi1_matrix + u*e_tilde
  res4P_ubar <- reverse.MHLT(psi.VARG,u1,u2,H,psi.parameterization=model)
  u1 <- gamma_matrix
  res4P_lbar <- reverse.MHLT(psi.VARG,u1,u2,H,psi.parameterization=model)

  xi1_3Dmatrix <- array(xi1,c(n.w,n.e,H))
  A_P_ubar <- - xi1_3Dmatrix + res4P_ubar$A
  A_P_lbar <- - xi1_3Dmatrix + res4P_lbar$A

  xi0_3Dmatrix  <- (- xi0) * array((1:H)%x%rep(1,n.e),c(1,n.e,H))
  B_P_ubar <- xi0_3Dmatrix + res4P_ubar$B
  B_P_lbar <- xi0_3Dmatrix + res4P_lbar$B

  return(list(A_P_ubar=A_P_ubar,
              A_P_lbar=A_P_lbar,
              B_P_ubar=B_P_ubar,
              B_P_lbar=B_P_lbar))
}


compute_proba_def_VARG <- function(model,H=10,indic_delta1 = NaN,
                                   X=NaN){
  # model contains, in particular, the parameterization of the VARG model.
  # If indic_delta1 is an integer, it indicates the component of w_t that
  #  corresponds to delta_{e=1}. If NaN, delta_{e=1} is supposed to correspond
  #  to the first non-zero entry of vector nu.

  if(is.na(indic_delta1[1])){
    indic_delta1 <- which(nu==0)[1]
  }
  n.w <- dim(model$beta)[1]
  n.e <- n.w - indic_delta1 + 1
  n.y <- n.w - n.e

  e_tilde <- rbind(matrix(0,n.y,n.e),
                   diag(n.e))

  u <- -10000

  u1 <- u*e_tilde
  u2 <- u*e_tilde
  res <- reverse.MHLT(psi.VARG,u1,u2,H,psi.parameterization=model)

  A <- res$A
  B <- res$B

  if(is.na(X[1])){
    PD <- NaN
  }else{
    T <- dim(X)[1]
    vec1T <- matrix(1,T,1)
    PD <- array(NaN,c(T,n.e,H))
    for(h in 1:H){
      PD[,,h] <- 1 - exp(X %*% A[,,h] + vec1T %*% B[1,,h])
    }
  }

  return(list(A = A, B = B,
              PD = PD))
}


compute_uncondmean_VARG <- function(model){
  alpha <- model$alpha
  nu    <- model$nu
  mu    <- model$mu
  beta  <- model$beta

  n <- dim(beta)[1]
  mu_matrix <- matrix(mu,n,n)
  rho <- mu_matrix * beta

  E <- solve(diag(n) - rho) %*% (mu * (alpha + nu))

  return(E)
}

simul_VARG <- function(model,nb_periods,w0=NaN){
  alpha <- model$alpha
  nu    <- model$nu
  mu    <- model$mu
  beta  <- model$beta

  n <- dim(beta)[1] # dimension of state vector

  if(is.na(w0[1])){
    w0 <- compute_uncondmean_VARG(model)
  }

  w     <- w0
  all_w <- w
  for(t in 1:nb_periods){
    param.pois <- rpois(n, alpha + beta %*% w)
    w <- rgamma(n,shape = param.pois + nu, scale = mu)
    all_w <- cbind(all_w,w)
  }

  return(all_w)
}

make_VARG_Q <- function(model){
  modelQ <- model
  modelQ$mu    <- model$mu / (1 - model$alpha_w * model$mu)
  modelQ$alpha <- model$alpha * modelQ$mu/model$mu
  modelQ$beta  <- diag(c(modelQ$mu/model$mu)) %*% model$beta
  return(modelQ)
}

prices_CDS_RFV_VARG <- function(model,H=10,indic_delta1 = NaN,
                                X=NaN){
  # X is a matrix of dimension T x n.w of state vector values.
  # Prices are computed directly for all defaultable entities considered
  #      in the VARG model (number = n.e).
  # If indic_delta1 is an integer, it indicates the component of w_t that
  #      corresponds to delta_{e=1}. If NaN, delta_{e=1} is supposed to
  #      correspond to the first non-zero entry of vector nu.
  # Currency: This function can handle the pricing of CDS written in a currency
  #           that is different from the currency of denomination of the
  #           underyling bond. For that model$mu_s0 and model$mu_s1 have
  #           to be different from 0.

  if(is.na(indic_delta1[1])){
    indic_delta1 <- which(nu==0)[1]
  }
  n.w <- dim(model$beta)[1]
  n.e <- n.w - indic_delta1 + 1
  n.y <- n.w - n.e

  mu_R0 <- model$mu_R0
  mu_R1 <- matrix(model$mu_R1,ncol=1)

  if(is.na(model$mu_s0)){
    mu_s0 <- 0
    mu_s1 <- matrix(0,n.w,1)
  }else{
    mu_s0 <- model$mu_s0
    mu_s1 <- matrix(model$mu_s1,ncol=1)
  }

  res_P_bar_0 <- compute_P_bar_VARG(model,
                                    mu_s1,
                                    H=H,
                                    indic_delta1 = indic_delta1)
  res_P_bar_muR1 <- compute_P_bar_VARG(model,
                                       mu_s1 - mu_R1,
                                       H=H,
                                       indic_delta1 = indic_delta1)

  if(is.na(X[1,1])){
    CDS_spreads <- NaN
  }else{
    T <- dim(X)[1]
    vec1T <- matrix(1,T,1)

    CDS_spreads <- array(NaN,c(T,n.e,H))

    numerator <- 0
    denominat <- 0

    matrix_mu_R0 <- t(matrix(mu_R0,n.e,T))

    for(h in 1:H){
      P_lb_0 <- exp(vec1T %*% res_P_bar_0$B_P_lbar[1,,h] +
                      X %*% res_P_bar_0$A_P_lbar[,,h])
      P_ub_0 <- exp(vec1T %*% res_P_bar_0$B_P_ubar[1,,h] +
                      X %*% res_P_bar_0$A_P_ubar[,,h])
      P_lb_muR1 <- exp(vec1T %*% res_P_bar_muR1$B_P_lbar[1,,h] +
                         X %*% res_P_bar_muR1$A_P_lbar[,,h])
      P_ub_muR1 <- exp(vec1T %*% res_P_bar_muR1$B_P_ubar[1,,h] +
                         X %*% res_P_bar_muR1$A_P_ubar[,,h])

      numerator <- numerator + exp(mu_s0) * (P_lb_0 - P_ub_0 -
        exp(-matrix_mu_R0) * (P_lb_muR1 - P_ub_muR1))
      denominat <- denominat + exp(mu_s0) * P_ub_0

      CDS_spreads[,,h] <- numerator/denominat
    }
  }

  return(list(
    res_P_bar_0    = res_P_bar_0,
    res_P_bar_muR1 = res_P_bar_muR1,
    CDS_spreads    = CDS_spreads
  ))
}


