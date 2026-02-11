
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
  A.h_1 <- 0
  B.h_1 <- 0
  n <- dim(u1)[1]
  k <- dim(u1)[2]
  A <- array(NaN,c(n,k,H))
  B <- array(NaN,c(1,k,H))
  for(i in 1:H){
    if(i==1){u <- u1}else{u <- u2}
    psi.u <- psi(u + A.h_1,psi.parameterization)
    A.h <- matrix(psi.u$a,n,k)
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
  if(dim(u)[1]==1){
    aux <- matrix(apply(u,2,function(x){x %x% x}),ncol=1)
  }else{
    aux <- t(apply(u,2,function(x){x %x% x}))
  }
  b <- t(u) %*% mu + .5 * aux %*% c(Sigma)
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




simul.ARG_OLD <- function(nb.sim,mu,nu,rho,alpha=0,w0=NaN){
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

simul.ARG <- function(nb.sim,mu,nu,rho,alpha=0,w0=NaN,nb.replic=1){
  # This function simulates an ARG or an ARG0 process
  # If nb.replic>1, then we simulate several paths in parallel

  unc.mean <- (alpha + nu) * mu/(1-rho)
  # set starting point:
  if(is.na(w0)){
    w_0 <- unc.mean
  }else{
    w_0 <- w0
  }

  W <- matrix(NaN,nb.sim,nb.replic)
  W[1,] <- w_0
  w <- w_0
  for(t in 2:nb.sim){
    z <- rpois(nb.replic,rho*w/mu + alpha)
    w <- mu * rgamma(nb.replic,shape=nu+z)
    W[t,] <- w
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
  # psi admits 2 arguments: model and u.
  #    (model is a list of parameters; u is the vector of points
  #     where one wants to evaluate the pdf.)

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
  # Impose appropriate format (column vector for x, row vector for gamma):
  x     <- matrix(x,ncol=1)
  gamma <- matrix(gamma,nrow=1)

  u <- c(1i*x)
  psi_eval <- matrix(psi(u,model),length(x),length(gamma))
  dx <- matrix(x-c(0,x[1:length(x)-1]),length(x),length(gamma))
  fx <- Im(psi_eval*exp(-1i*x %*% gamma))/matrix(x,length(x),length(gamma))*dx

  f  <- 1/2 - 1/pi * apply(fx,2,sum)

  f.cdf <- pmax(pmin(f,1),0)
  for(i in 2:length(f)){
    if(f.cdf[i]<f.cdf[i-1]){
      f.cdf[i] <- f.cdf[i-1]}}

  return(f.cdf)
}

# simul.VAR <- function (Model, nb.sim, x0 = NaN)
# {
#   n <- dim(Model$Phi)[1]
#   if (is.na(x0[1])) {
#     x0 <- solve(diag(n) - Model$Phi) %*% Model$mu
#   }
#   X <- c(x0)
#   x <- x0
#   for (t in 2:nb.sim) {
#     x <- Model$mu + Model$Phi %*% x + Model$Sigma %*% rnorm(n)
#     X <- rbind(X, c(x))
#   }
#   return(X)
# }

compute_expect_variance <- function(psi,model,du = 1e-06){
  # This function computes the conditional and unconditional expectations and
  # varances of an affine process.
  # ----------------------------------------------------------------------------
  # The conditional mean is given by:
  # E_t(w_{t+1}) = mu + Phi.w_t
  # The conditional variance is given by:
  # vec[Var_t(w_{t+1})] = Theta0 + Theta1.w_t
  # ----------------------------------------------------------------------------

  if(is.null(model$n_w)){
    print("The 'model' input should contain 'n_w', the dimension of the state vector.")
    return(0)
  }else{
    n_w  <- model$n_w
  }

  # Computation of Ew, using the Laplace Transform (Ew = dPsi(u)/du at u=0) ----
  mu <- matrix(0,n_w,1)
  Phi <- matrix(0,n_w,n_w)
  Theta0 <- matrix(0,n_w*n_w,1)
  Theta1 <- matrix(0,n_w*n_w,n_w)

  for(i in 1:n_w){
    u_plus_i <- matrix(0,model$n_w,1)
    u_minu_i <- matrix(0,model$n_w,1)
    u_plus_i[i] <- du
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
      Theta0[n_w*(i-1)+j,1] <-
        (psi_plus_i_plus_j$b - psi_plus_i_minu_j$b - psi_minu_i_plus_j$b + psi_minu_i_minu_j$b)/
        (4*du^2)
      Theta1[n_w*(i-1)+j,]  <-
        (psi_plus_i_plus_j$a - psi_plus_i_minu_j$a - psi_minu_i_plus_j$a + psi_minu_i_minu_j$a)/
        (4*du^2)
      # print(cbind(psi_plus_i_plus_j$a,
      #             psi_plus_i_minu_j$a,
      #             psi_minu_i_plus_j$a,
      #             psi_minu_i_minu_j$a))
      # print(psi_plus_i_plus_j$a - psi_plus_i_minu_j$a - psi_minu_i_plus_j$a + psi_minu_i_minu_j$a)
      # return(0)
    }
  }

  # Compute unconditional expectation:
  Ew <- solve(diag(n_w) - Phi) %*% mu
  # Compute unconditional variance:
  Vw <- matrix(solve(diag(n_w^2) - Phi %x% Phi) %*% (Theta0 + Theta1 %*% Ew),
               n_w,n_w)

  return(list(
    mu = mu,
    Phi = Phi,
    Theta0 = Theta0,
    Theta1 = Theta1,
    Ew = Ew,
    Vw = Vw
  ))
}


compute_expect_variance_H <- function(VAR_representation,
                                      theta1,theta2=theta1,H){
  # The model is:
  # w_{t+1} = mu + Phi*w_{t} + eps_{t+1},
  # with E_t(eps_{t+1}) = 0 and vec[Var_t(eps_{t+1})]=Theta0 + Theta1*w_{t}.
  # In that case, we have
  # E_t(theta2*w_{t+1} + ... + theta2*w_{t+h-1} + theta1*w_{t+h}) =
  #      C0(theta2,theta1,h) + C1(theta2,theta1,h)*w_{t} and
  # vec[Var_t(theta2*w_{t+1} + ... + theta2*w_{t+h-1} + theta1*w_{t+h})] =
  #      Gamma0(theta2,theta1,h) + Gamma1(theta2,theta1,h)*w_{t}.

  mu  <- VAR_representation$mu
  Phi <- VAR_representation$Phi
  Theta0 <- VAR_representation$Theta0
  Theta1 <- VAR_representation$Theta1

  k <- dim(theta1)[1]
  n <- dim(theta1)[2]

  all_C0 <- array(NaN,c(k,1,H))
  all_C1 <- array(NaN,c(k,n,H))
  all_Gamma0 <- array(NaN,c(k^2,1,H))
  all_Gamma1 <- array(NaN,c(k^2,n,H))

  for(h in 1:H){
    if(h==1){
      C0_h <- theta1 %*% mu
      C1_h <- theta1 %*% Phi
      Gamma0_h <- (theta1 %x% theta1) %*% Theta0
      Gamma1_h <- (theta1 %x% theta1) %*% Theta1
    }else{
      C0_h <- (theta2 + C1_h_1) %*% mu + C0_h_1
      C1_h <- (theta2 + C1_h_1) %*% Phi

      aux <- (theta2 + C1_h_1) %x% (theta2 + C1_h_1)
      Gamma0_h <- Gamma0_h_1 + Gamma1_h_1 %*% mu + aux %*% Theta0
      Gamma1_h <- Gamma1_h_1 %*% Phi + aux %*% Theta1
    }
    C0_h_1 <- C0_h
    C1_h_1 <- C1_h
    Gamma0_h_1 <- Gamma0_h
    Gamma1_h_1 <- Gamma1_h

    all_C0[,,h] <- C0_h
    all_C1[,,h] <- C1_h
    all_Gamma0[,,h] <- Gamma0_h
    all_Gamma1[,,h] <- Gamma1_h
  }

  return(list(all_C0 = all_C0,
              all_C1 = all_C1,
              all_Gamma0 = all_Gamma0,
              all_Gamma1 = all_Gamma1))
}

# # Check
# modelVAR <- list(mu = matrix(c(1,0),2,1),
#                  Phi = matrix(c(.6,.2,-.3,.8),2,2),
#                  Sigma = diag(2),
#                  n_w=2)
#
# VAR_representation <- compute_expect_variance(psi.GaussianVAR,modelVAR)
# solve(diag(4) - modelVAR$Phi %x% modelVAR$Phi) %*% c(modelVAR$Sigma)
#
# theta1 <- matrix(rnorm(8),4,2)
# H <- 5
# res <- compute_expect_variance_H(VAR_representation,
#                                  theta1=theta1,H=H)


psi.VARG <- function(u,psi.parameterization){
  # Laplace transform of a VARG with conditionally independent components.
  # Conditionally on w_t,w_{t-1},..., the distribution of the jth component is:
  #  w_{j,t} ~ gamma_{nu_j}(alpha_j + beta_j'w_{t-1},mu_j),
  # In psi.parameterization, the vectors alpha, nu, mu gather the n values
  #  of the associated parameters. The rows of matrix beta are the beta_j'
  # If w_t is n-dimensional (of dimension n, say), u is of dimension n x k
  #    (i.e., we can compute k LT in parallel)
  alpha <- psi.parameterization$alpha
  beta  <- psi.parameterization$beta
  nu    <- psi.parameterization$nu
  mu    <- psi.parameterization$mu
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

compute_uncondvar_VARG <- function(model){
  alpha <- model$alpha
  nu    <- model$nu
  mu    <- model$mu
  beta  <- model$beta

  n <- dim(beta)[1]
  mu_matrix <- matrix(mu,n,n)
  rho <- mu_matrix * beta

  E       <- solve(diag(n) - rho) %*% (mu * (alpha + nu))
  Sigma_E <- diag(c(mu^2*(nu + 2*alpha) + 2*(mu_matrix^2*beta) %*% E))
  V       <- matrix(solve(diag(n^2) - rho %x% rho) %*% c(Sigma_E),n,n)

  return(V)
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
  w <- w0
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
  mu_R1 <- matrix(model$mu_R1,n.w,n.e)

  if(is.null(model$mu_s0)){
    mu_s0 <- 0
    mu_s1 <- matrix(0,n.w,n.e)
  }else{
    mu_s0 <- model$mu_s0
    mu_s1 <- matrix(model$mu_s1,n.w,n.e)
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
    matrix_mu_s0 <- t(matrix(mu_s0,n.e,T))

    for(h in 1:H){
      P_lb_0 <- exp(vec1T %*% res_P_bar_0$B_P_lbar[1,,h] +
                      X %*% res_P_bar_0$A_P_lbar[,,h])
      P_ub_0 <- exp(vec1T %*% res_P_bar_0$B_P_ubar[1,,h] +
                      X %*% res_P_bar_0$A_P_ubar[,,h])
      P_lb_muR1 <- exp(vec1T %*% res_P_bar_muR1$B_P_lbar[1,,h] +
                         X %*% res_P_bar_muR1$A_P_lbar[,,h])
      P_ub_muR1 <- exp(vec1T %*% res_P_bar_muR1$B_P_ubar[1,,h] +
                         X %*% res_P_bar_muR1$A_P_ubar[,,h])

      numerator <- numerator + exp(matrix_mu_s0) * (P_lb_0 - P_ub_0 -
                                                      exp(-matrix_mu_R0) * (P_lb_muR1 - P_ub_muR1))
      denominat <- denominat + exp(matrix_mu_s0) * P_ub_0

      CDS_spreads[,,h] <- numerator/denominat
    }
  }

  return(list(
    res_P_bar_0    = res_P_bar_0,
    res_P_bar_muR1 = res_P_bar_muR1,
    CDS_spreads    = CDS_spreads
  ))
}


# Top-down approach ============================================================

psi.w.TopDown <- function(u,model,psi.y=psi.y.TopDown){
  # This function evaluates the LT of the state vector w_t in the context of the
  # top-down approach of Gourieroux, Monfort, Mouabbi, and Renne.
  # "model" is a list containing the model parameters.
  # "psi.y" is the LT of y_{t+1} conditionally on w_t.

  alpha.n <- matrix(model$alpha.n,ncol=1)
  beta.ny <- model$beta.ny
  beta.nn <- model$beta.nn

  J   <- dim(alpha.n)[1] # number of segments in the economy.
  n.w <- dim(u)[1]
  n.y <- n.w - J

  u.y <- u[1:n.y,]
  u.n <- u[(n.y+1):(n.y+J),]

  ab_y <- psi.y.TopDown(u.y + t(beta.ny) %*% (exp(u.n) - 1),model)

  a.wy <- ab_y$a.yy
  a.wn <- ab_y$a.yn + t(beta.nn) %*% (exp(u.n) - 1)
  b.w  <- ab_y$b    + t(exp(u.n) - 1) %*% alpha.n

  return(list(b    = b.w,
              a.wy = a.wy,
              a.wn = a.wn,
              a    = rbind(a.wy,a.wn)))
}

psiQ.w.TopDown <- function(u,model,psi.y=psi.y.TopDown){
  # Same as psiQ.w.TopDown, but under the risk-neutral (Q) measure

  k <- dim(u)[2]
  u_adjusted <- cbind(u + model$alpha.m %*% matrix(1,1,k),model$alpha.m)

  res <- psi.w.TopDown(u_adjusted,model,psi.y)
  b.w  <- res$b
  a.wy <- res$a.wy
  a.wn <- res$a.wn

  a.wy <- a.wy[,1:k] - a.wy[,k+1] %*% matrix(1,1,k)
  a.wn <- a.wn[,1:k] - a.wn[,k+1] %*% matrix(1,1,k)
  b.w  <- b.w[1:k,]  - b.w[k+1]

  return(list(b    = b.w,
              a.wy = a.wy,
              a.wn = a.wn,
              a    = rbind(a.wy,a.wn)))
}

psi.y.TopDown <- function(u.y,model){
  # This function returns the LT of y_{t+1}|w_t.
  # It can be evaluated for different u.y's, i.e., u.y can be a matrix.
  # In this example, the conditional distribution of y_t is a compound
  #   Poisson-gamma (as in ARG processes).

  beta.yy <- model$beta.yy
  beta.yn <- model$beta.yn

  alpha.y <- model$alpha.y
  nu.y    <- model$nu.y
  mu.y    <- model$mu.y
  n.y <- dim(u.y)[1]
  k   <- dim(u.y)[2]
  mu_matrix <- matrix(mu.y,n.y,k)
  a.yy <- t(beta.yy) %*% ((u.y * mu_matrix)/(1 - u.y * mu_matrix))
  a.yn <- t(beta.yn) %*% ((u.y * mu_matrix)/(1 - u.y * mu_matrix))
  b.y <- t((u.y * mu_matrix)/(1 - u.y * mu_matrix)) %*% alpha.y -
    t(log(1 - u.y * mu_matrix)) %*% nu.y
  return(list(a.yy = a.yy,
              a.yn = a.yn,
              b    = b.y))
}

simul.TopDown <- function(model,nb_periods,W0=NaN){
  # This function simulates the Top Down model. The conditional distribuion of
  # y_t is compound Poisson-Gamma (as in the ARG).
  beta.yy <- as.matrix(model$beta.yy)
  beta.yn <- as.matrix(model$beta.yn)
  beta.nn <- as.matrix(model$beta.nn)
  beta.ny <- as.matrix(model$beta.ny)

  alpha.n <- matrix(model$alpha.n,ncol=1)

  alpha.y <- matrix(model$alpha.y,ncol=1)
  nu.y    <- matrix(model$nu.y,ncol=1)
  mu.y    <- matrix(model$mu.y,ncol=1)

  n.y <- dim(alpha.y)[1]
  J   <- dim(beta.nn)[1]
  n.w <- n.y + J

  if(is.na(W0)[1]){
    EV <- compute_expect_variance(psi.w.TopDown,model)
    Ew <- EV$Ew
    W0 <- matrix(Ew,ncol=1)
  }
  W <- matrix(W0,ncol=1)
  y <- matrix(W0[1:n.y],ncol=1)
  n <- matrix(W0[(n.y+1):n.w],ncol=1)
  all_y <- y
  all_n <- n

  for(t in 2:nb_periods){
    # y part:
    param.pois <- rpois(n.y, alpha.y + beta.yy %*% y + beta.yn %*% n)
    y <- rgamma(n.y,shape = param.pois + nu.y, scale = mu.y)
    all_y <- cbind(all_y,y)
    # n part:
    n <- rpois(J, alpha.n + beta.ny %*% y + beta.nn %*% n)
    all_n <- cbind(all_n,n)
  }

  return(list(
    EV = EV,
    all_y = all_y,
    all_n = all_n,
    all_w = rbind(all_y,
                  all_n)
  ))
}


compute_H_bar_TopDown <- function(model,gamma,H=10,W=NaN,
                                  indic_Q = TRUE){
  # gamma has to be a vector with the same dimension as w_t
  # model contains the parameterization of the VARG model.
  # Prices are computed for all segments.

  gamma <- matrix(gamma,ncol=1)
  xi1 <- matrix(model$xi1,ncol=1)
  xi0 <- model$xi0

  epsilon <- 10^(-7) # used to numerically compute gradient

  # Build u1 and u2 -------------
  # (1st column: for H_upper; 2nf column: for H_lower,
  #  3rd and 4th columns: are the same, for H0)
  u1 <- cbind(gamma*epsilon,
              gamma*0,
              gamma*0,
              gamma*0)
  u2 <- cbind(gamma*epsilon - xi1,
              gamma*epsilon - xi1,
              gamma*0       - xi1,
              gamma*0       - xi1)

  if(indic_Q){
    psi <- psiQ.w.TopDown
  }else{
    psi <- psi.w.TopDown
  }
  varphi <- reverse.MHLT(psi,
                         u1 = u1,
                         u2 = u2,
                         H = H,model)
  nw <- dim(varphi$A)[1]

  da_du <- (varphi$A[,1:2,] - varphi$A[,3:4,])/epsilon
  db_du <- (varphi$B[,1:2,] - varphi$B[,3:4,])/epsilon
  db_du <- array(db_du,c(1,2,H))

  a_exp <- varphi$A[,1:2,] - array(xi1,c(nw,2,H))
  b_exp <- array(varphi$B[,1:2,],c(1,2,H)) -
    array(xi0*(1:H)%x%c(1,1),c(1,2,H))

  # Risk-free bond prices:
  a_exp_rf <- varphi$A[,3:4,] - array(xi1,c(nw,2,H))
  b_exp_rf <- array(varphi$B[,3:4,],c(1,2,H)) -
    array(xi0*(1:H)%x%c(1,1),c(1,2,H))
  a_exp_rf <- a_exp_rf[,1,]
  b_exp_rf <- b_exp_rf[,1,]

  if(!is.na(W[1])){
    if(dim(W)[2]!=nw){
      # to make sure W is T x nw
      W <- t(W)
    }
    T <- dim(W)[1]
    H_upper <- matrix(0,T,H)
    H_lower <- matrix(0,T,H)
    B       <- matrix(0,T,H)
    for(h in 1:H){
      H_upper[,h] <- exp(b_exp[1,1,h] + W %*% a_exp[,1,h]) *
        (db_du[1,1,h] + W %*% da_du[,1,h])
      H_lower[,h] <- exp(b_exp[1,2,h] + W %*% a_exp[,2,h]) *
        (db_du[1,2,h] + W %*% da_du[,2,h])
      B[,h] <- exp(b_exp_rf[h] + W %*% a_exp_rf[,h])
    }
  }else{
    H_upper <- NaN
    H_lower <- NaN
    B       <- NaN
  }
  return(list(da_du = da_du,
              db_du = db_du,
              a_exp = a_exp,
              b_exp = b_exp,
              a_exp_rf = a_exp_rf,
              b_exp_rf = b_exp_rf,
              H_upper = H_upper,
              H_lower = H_lower,
              B = B))
}


compute_CDS_TopDown <- function(model,
                                H=10,
                                W,
                                RR=.5,
                                indic_Q=TRUE){

  alpha.n <- matrix(model$alpha.n,ncol=1)

  J   <- length(model$alpha.n) # number of segments in the economy.
  n.y <- length(model$alpha.y) # number of y variables.
  n.w <- J + n.y

  if(dim(W)[2]!=n.w){# to make sure W is of dimension T x n.w
    W <- t(W)
  }
  T <- dim(W)[1]

  CDS <- array(0,c(T,J,H))

  for(j in 1:J){
    I_j <- model$I[j]
    ej_tilde <- matrix(0,n.w,1)
    ej_tilde[n.y+j] <- 1
    res <- compute_H_bar_TopDown(model,ej_tilde,H=H,W,
                                 indic_Q = indic_Q)
    numerator <- 0
    denominat <- 0
    for(h in 1:H){
      numerator <- numerator + res$H_upper[,h] - res$H_lower[,h]
      denominat <- denominat + I_j * res$B[,h] - res$H_lower[,h]
      CDS[,j,h] <- (1-RR)*(numerator/denominat)
    }
  }

  return(CDS)
}


varphi4G_TopDown <- function(x,parameterization,H){

  model   <- parameterization$model
  gamma   <- parameterization$gamma
  indic_Q     <- parameterization$indic_Q # if FALSE, use physical LT, otherwise Q
  indic_upper <- parameterization$indic_upper # if TRUE, compute G_upper, otherwise G_lower

  v <- parameterization$v # determine thresholds (dimension nw x 1)
  i.v.x <- matrix(v,ncol=1) %*% matrix(c(1i*x),nrow=1)

  nw <- dim(i.v.x)[1]
  q  <- dim(i.v.x)[2]

  u2 <- matrix(gamma - model$xi1, nw, q) + i.v.x
  u1 <- i.v.x
  if(indic_upper){
    u1 <- u1 + matrix(gamma, nw, q)
  }

  if(indic_Q){
    psi <- psiQ.w.TopDown
  }else{
    psi <- psi.w.TopDown
  }

  res_reverse <- reverse.MHLT(psi,
                              u1 = u1,
                              u2 = u2,
                              H = H,
                              psi.parameterization = model)

  A <- res_reverse$A
  A <- A - array(matrix(1,q*H,1) %x% matrix(model$xi1,ncol=1),c(nw,q,H))
  B <- res_reverse$B
  B <- B - model$xi0*array(matrix((1:H),ncol=1) %x% matrix(1,q,1),c(1,q,H))

  return(list(A = res_reverse$A,
              B = matrix(res_reverse$B,nrow=1)))
}


truncated.payoff <- function(W, # values of w_t
                             b.matrix, # thresholds (H x k, 1 row per maturity)
                             varphi,
                             parameterization,
                             max_x = 500,
                             dx_statio = .1,
                             min_dx = 1e-05,
                             nb_x1 = 1000){
  # This function computes:
  #       E[exp(u0·w_t+...+uh·w_{t+h})·1_{v0·w_t+...+vh·w_{t+h} < b}],
  # for different horizons and values of b.
  # b.matrix is organized as follows: H * k; the h^th row gives the b values for
  #       horizon h.
  # "W" is of dimension T x n.w; it contains values of the state vector w
  #       at which we want to evaluate the prices.
  # "varphi" is a function that takes three arguments:
  #       (i) u and
  #       (ii) parameterization (that contains at least "model"),
  #       (iii) H (max maturity);
  # it returns matrices A(u) and B(u) that are such that the "price" of
  #       the payoff is given by exp(A(u)'w + B(u)),
  #       where A(u) is of dimension n.u x n.w,
  #       and B(u) is of dimension n.u x 1
  # NOTE: It is varphi that determines the ui's and vi's.

  if(is.null(parameterization$model$n_w)){
    print("Error: 'parameterization' should include 'model', itself including 'n_w' (size of state vector).")
    return(1)
  }

  n_w <- parameterization$model$n_w

  dx1 <- exp(seq(log(min_dx),log(dx_statio),length.out=nb_x1))
  max_x1 <- tail(cumsum(dx1),1)
  nb_x2 <- (max_x - max_x1)/dx_statio
  dx <- c(dx1,rep(dx_statio,nb_x2))
  x <- matrix(cumsum(dx),ncol=1)

  nb_x <- length(x)

  H <- dim(b.matrix)[1]

  res_varphi  <- varphi(x,parameterization,H)
  res_varphi0 <- varphi(0,parameterization,H)

  # Adjust for current short-term interest rate:
  B <- matrix(res_varphi$B,1,nb_x*H)
  A <- matrix(res_varphi$A,n_w,nb_x*H)
  B0 <- matrix(res_varphi0$B,nrow=1)
  A0 <- matrix(res_varphi0$A,n_w,H)

  if(dim(W)[2]!=n_w){# to make sure W is of dimension T x n.w
    W <- t(W)
  }
  T <- dim(W)[1]

  psi_eval  <- exp(t(A)  %*% t(W) + matrix(B, nb_x*H,T))
  psi_eval0 <- Re(exp(t(A0) %*% t(W) + matrix(B0,H,T)))

  # modify format to have them nb_x x (T*H):
  Psi_eval  <- matrix(NaN,nb_x,T*H)
  for(h in 1:H){
    Psi_eval[,(T*(h-1)+1):(T*h)]  <- psi_eval[(nb_x*(h-1)+1):(nb_x*h),]
  }

  dx <- matrix(x - c(0,x[1:(nb_x-1)]),nb_x,T*H)
  x_mid <- matrix(x,nb_x,T*H)

  all_prices <- matrix(0,T*H,dim(b.matrix)[2])

  count_b <- 0
  for(b in 1:dim(b.matrix)[2]){
    count_b <- count_b + 1
    aux <- Im(Psi_eval *
                exp(-1i*matrix(x,ncol=1) %*%
                      (matrix(b.matrix[,count_b],nrow=1) %x% matrix(1,1,T))))
    fx <- aux/x_mid * dx
    f  <- 1/2*c(t(psi_eval0)) - 1/pi * apply(fx,2,sum)
    all_prices[,count_b] <- f
  }

  All_prices <- array(NaN,c(T,dim(b.matrix)[2],H))
  for(h in 1:H){
    All_prices[,,h] <- all_prices[(T*(h-1)+1):(T*h),]
  }

  return(All_prices)
}


compute_G <- function(model,W,
                      gamma,v,
                      b.matrix, # attachment/detachment points
                      H, # maturity
                      indic_Q = TRUE, # computed under Q.
                      indic_upper = TRUE, # otherwise, computed G_lower
                      max_x = 2000,
                      dx_statio = 2,
                      min_dx = 10^(-6),
                      nb_x1=1000){
  # This computes functions G_lower and G_upper

  parameterization <- list(model=model,
                           indic_Q = indic_Q,
                           indic_upper = indic_upper,
                           gamma = gamma,
                           v = v)

  # P <- truncated.payoff(W,
  #                       v,vector.of.b,
  #                       H,
  #                       varphi = varphi4G_TopDown,
  #                       parameterization,
  #                       max_x = max_x,
  #                       dx_statio = dx_statio,
  #                       min_dx = min_dx,
  #                       nb_x1=nb_x1)
  P <- truncated.payoff(W,
                        b.matrix,
                        varphi = varphi4G_TopDown,
                        parameterization,
                        max_x = max_x,
                        dx_statio = dx_statio,
                        min_dx = min_dx,
                        nb_x1=nb_x1)

  return(P)
}

compute_tranche <- function(model,W,
                            j,H,
                            vector.of.a, # attachment/detachment points
                            RR = 0.5, # recovery rate
                            max_x = 500,
                            dx_statio = .1,
                            min_dx = 1e-05,
                            nb_x1 = 1000,
                            indic_Q = TRUE # if FALSE; then everyting computed under P
){
  # This function computes the price of credit tranche products in
  #     the context of the top-down credit approach.
  # It does that for one segment only (segment j).
  # vector.of.a contains the attachment/detachment points. For instance,
  #     if vector.of.a =[0,.03,.06], then we consider 2 tranches:
  #     [0,3%] and [3%,6%].
  # W is a matrix of dimension T x n.w containing the consider values
  #     of the state vector.

  n_w <- model$n_w
  if(dim(W)[2]!=n_w){# to make sure W is of dimension T x n.w
    W <- t(W)
  }
  T <- dim(W)[1]

  vector.of.a.bar <- vector.of.a * model$I[j] / (1 - RR)
  nb_a <- length(vector.of.a)

  # Prepare arrays of a.bar and b.bar
  # b.bar:
  aux <- vector.of.a.bar[2:nb_a]
  aux <- matrix(1,T,1) %*% matrix(aux,1,nb_a-1)
  array.b.bar <- array(aux,c(T,nb_a-1,H))
  # a.bar:
  aux <- vector.of.a.bar[1:(nb_a-1)]
  aux <- matrix(1,T,1) %*% matrix(aux,1,nb_a-1)
  array.a.bar <- array(aux,c(T,nb_a-1,H))

  n.w <- model$n_w
  J <- dim(as.matrix(model$beta.nn))[1]
  n.y <- n.w - J

  ej_tilde <- matrix(0,n.w,1)
  ej_tilde[n.y + j] <- 1

  G_upper_0_ej <- compute_G(model,
                            W=W,
                            gamma = 0*ej_tilde,
                            v = ej_tilde,
                            b.matrix=t(matrix(vector.of.a.bar,
                                              length(vector.of.a.bar),H)), # attachment/detachment points
                            H=H,
                            indic_Q = indic_Q,
                            indic_upper = TRUE,
                            max_x = max_x,dx_statio = dx_statio,
                            min_dx = min_dx,nb_x1 = nb_x1)
  G_lower_0_ej <- compute_G(model,
                            W=W,
                            gamma = 0*ej_tilde,
                            v = ej_tilde,
                            b.matrix=t(matrix(vector.of.a.bar,
                                              length(vector.of.a.bar),H)), # attachment/detachment points
                            H=H,
                            indic_Q = indic_Q,
                            indic_upper = FALSE,
                            max_x = max_x,dx_statio = dx_statio,
                            min_dx = min_dx,nb_x1 = nb_x1)

  u <- 10^(-4)
  G_upper_uej_ej <- compute_G(model,
                              W=W,
                              gamma = u*ej_tilde,
                              v = ej_tilde,
                              b.matrix=t(matrix(vector.of.a.bar,
                                                length(vector.of.a.bar),H)), # attachment/detachment points
                              H=H,
                              indic_Q = indic_Q,
                              indic_upper = TRUE,
                              max_x = max_x,dx_statio = dx_statio,
                              min_dx = min_dx,nb_x1 = nb_x1)
  G_lower_uej_ej <- compute_G(model,
                              W=W,
                              gamma = u*ej_tilde,
                              v = ej_tilde,
                              b.matrix=t(matrix(vector.of.a.bar,
                                                length(vector.of.a.bar),H)), # attachment/detachment points
                              H=H,
                              indic_Q = indic_Q,
                              indic_upper = FALSE,
                              max_x = max_x,dx_statio = dx_statio,
                              min_dx = min_dx,nb_x1 = nb_x1)

  dG_upper <- (G_upper_uej_ej - G_upper_0_ej)/u
  dG_lower <- (G_lower_uej_ej - G_lower_0_ej)/u

  G_upper_0_ej_a <- G_upper_0_ej[,1:(nb_a-1),]
  G_upper_0_ej_b <- G_upper_0_ej[,2:nb_a,]

  dG_upper_a <- dG_upper[,1:(nb_a-1),]
  dG_upper_b <- dG_upper[,2:nb_a,]

  dG_lower_a <- dG_lower[,1:(nb_a-1),]
  dG_lower_b <- dG_lower[,2:nb_a,]

  numerator <- dG_upper_b - dG_upper_a - (dG_lower_b - dG_lower_a)
  denominat <- dG_upper_a - array.a.bar*G_upper_0_ej_a +
    array.b.bar*G_upper_0_ej_b - dG_upper_b

  for(h in 2:H){
    numerator[,,h] <- numerator[,,h-1] + numerator[,,h]
    denominat[,,h] <- denominat[,,h-1] + denominat[,,h]
  }

  tranche_prices <- numerator/denominat

  return(list(tranche_prices=tranche_prices,
              numerator=numerator,
              denominat=denominat))
}


compute_affine_payoff_TopDown <- function(model,gamma,H,
                                          indic_Q = TRUE,
                                          W=NaN){
  # This functions computes the date-t price of exp(gamma'w_{t+h}) in
  # the context of the top-down approach.

  n_w <- model$n_w

  u2 <- matrix(- model$xi1, ncol=1)
  u1 <- matrix(gamma,n_w,1)

  if(indic_Q){
    psi <- psiQ.w.TopDown
  }else{
    psi <- psi.w.TopDown
  }

  res_reverse <- reverse.MHLT(psi,
                              u1 = u1,
                              u2 = u2,
                              H = H,
                              psi.parameterization = model)

  A <- matrix(res_reverse$A,n_w,H)
  B <- matrix(res_reverse$B,1,H)

  A <- A - matrix(model$xi1,ncol=1) %*% matrix(1,1,H)
  B <- B - matrix(model$xi0,ncol=1) %*% matrix(1:H,1,H)

  a <- A/t(matrix(1:H,H,n_w))
  b <- B/t(matrix(1:H,H,1))

  if(!is.na(W[1])){
    # make sure W of the size T x n_w
    if(dim(W)[1]==n_w){
      W <- t(W)
    }
    T <- dim(W)[1]
    P <- exp(W %*% A + matrix(1,T,1) %*% B)
    r <- -log(P)/t(matrix(1:H,H,T))
  }else{
    P <- NaN
    r <- NaN
  }

  return(list(A = A,
              B = B,
              a = a,
              b = b,
              P = P,
              r = r))
}

compute_bond_price_Ratings <- function(model,Y,H,psi){
  # This function computes bond prices in a context featuring
  #     credit ratings migration.
  # H is the maximum maturity considered.
  # Y is the matrix of y_t realizations. Dimension: T x n_y.
  # !!: model$kappa0 is of dimension 1 x K, where K is the number of
  #     ratings (including "Default"). model$kappa1 is of dimension n_y x K.
  #     The bottom entries of model$kappa0 and model$kappa1 are zeros.
  # "psi" is a function computing the Laplace transform of y_t. If this function
  #     needs arguments (on top of u), it has to be coded as a list
  #     named "psi.parameterization" included in "model".

  # Specification of stochastic migration probabilities:
  kappa0 <- matrix(model$kappa0,nrow=1) # of dimension 1 x K. Last entry should be 0.
  kappa1 <- model$kappa1 # of dimension n_y x K. Last row should be 0.
  K <- length(kappa0) # number of ratings
  n_y <- dim(kappa1)[1]

  # Specification of short-term rate:
  xi0 <- model$xi0
  xi1 <- matrix(model$xi1,ncol=1)

  u1 <- - kappa1
  u2 <- - kappa1 - xi1 %*% matrix(1,1,K)
  res_k <- reverse.MHLT(psi,u1=u1,u2=u2,H,model$psi.parameterization)

  # Adjust for xi0, kappa0, and the first xi1:
  B <- matrix(res_k$B,1,K*H) -
    matrix(1:H,nrow=1) %x% (matrix(1,1,K) * xi0) -
    matrix(1:H,nrow=1) %x% kappa0
  A <- matrix(res_k$A,n_y,K*H) - matrix(xi1,n_y,K*H)

  if(dim(Y)[2]!=n_y){# to make sure Y is of dimension T x n_y
    Y <- t(Y)
  }
  T <- dim(Y)[1]

  psi_eval <- exp(t(A)  %*% t(Y) + matrix(B, K*H,T))
  # organized this way: each line for a given date.
  # And then: h=1, k=1,...,K / h=2, k=1,...,K / h=3, etc.

  # Re-organize in array format (T x K x H):
  Psi_eval <- array(t(psi_eval),c(T,K,H))

  # Risk-free bond price:
  RF_bond_price <- Psi_eval[,K,] # of dimension T x H

  V_1 <- solve(model$V)
  V_1jK <- V_1[,K]
  V_ijV_1jK <- model$V * (matrix(1,K,1) %*% matrix(V_1jK,nrow=1)) # K x K

  Risky_bond_price   <- array(NaN,c(T,K-1,H))
  Risky_bond_yields  <- array(NaN,c(T,K-1,H))
  Risky_bond_spreads <- array(NaN,c(T,K-1,H))
  for(h in 1:H){
    bond_prices_h <- - V_ijV_1jK[1:(K-1),1:(K-1)] %*% t(Psi_eval[,1:(K-1),h])
    Risky_bond_price[,,h] <- t(bond_prices_h)
    Risky_bond_yields[,,h] <- - log(Risky_bond_price[,,h])/h
  }

  # compute yields:
  RF_bond_yields <- - log(RF_bond_price) * (matrix(1,T,1) %*% matrix(1/(1:H),nrow=1))

  # compute credit spreads:
  for(k in 1:(K-1)){
    Risky_bond_spreads[,k,] <- Risky_bond_yields[,k,] - RF_bond_yields
  }

  return(list(RF_bond_price = RF_bond_price,
              RF_bond_yields = RF_bond_yields,
              Risky_bond_price = Risky_bond_price,
              Risky_bond_yields = Risky_bond_yields,
              Risky_bond_spreads = Risky_bond_spreads))
}


simul.rating.migration <- function(model,Y,tau_ini=1){
  # This function simulates a trajectory of ratings for a given model and
  # a trajectory of the exogenous variable Y (of dimension T x n_y).

  T <- dim(Y)[1]

  tau_t <- tau_ini
  all_tau <- tau_t

  for(t in 2:T){
    P <-
      model$V %*%
      diag(c(exp(-t(model$kappa0) -  t(model$kappa1) %*% c(Y[t,])))) %*%
      solve(model$V)

    cumPi <- cumsum(P[tau_t,])
    tau_t <- which(runif(1)<cumPi)[1]

    all_tau <- c(all_tau,tau_t)
  }
  return(all_tau)
}


compute_AB_classical <- function(xi0,xi1,
                                 kappa0=NaN,kappa1=NaN,
                                 H, # maximum maturity
                                 psi,psi.parameterization){
  # This function employs recursive formulas to price risk-free and defaultable
  # bonds.
  # If the kappa's are NaN (default), the function prices
  #    credit-risk-free bonds only.
  # The risk-free short-term rate is r_t = xi0 + xi1'y_t.
  # The default intensity is lambda_t = kappa0 + kappa1'y_t.
  # kappa0 and kappa1 can correspond to different entities, that is, kappa0 can be
  #    of length E, an kappa1 of dimension n_y x E.
  # The function returns a list containing A, B, a, and b.
  # Bond prices are such that B_{t,h} = exp(B_h + A_h'y_t), where
  #    B_h is a scalar and A_h is a (n_y x 1) vector.
  # Bond yields are of the form b_h + a_h'y_t.

  n_y <- length(xi1)

  if(!is.na(kappa0[1])){
    # Pricing of risk-free & defaultable bonds
    E <- dim(kappa1)[2]
    u1 <- cbind(0,-kappa1)
    u2 <- cbind(-xi1,-matrix(xi1,n_y,E)-kappa1)
    xi0_kappa0 <- c(xi0,kappa0)
  }else{
    # Pricing of risk-free bonds only
    u1 <- matrix(rep(0,n_y),ncol=1)
    u2 <- matrix(-xi1,ncol=1)
    xi0_kappa0 <- xi0
  }

  res_reverse <- reverse.MHLT(psi = psi,
                              u1 = u1,
                              u2 = u2,
                              H = H,
                              psi.parameterization = psi.parameterization)
  A <- res_reverse$A
  B <- res_reverse$B
  # Adjust for xi0, kappa0, and the first xi1'y_t:

  a <- 0*A
  b <- 0*B

  for(h in 1:H){
    A[,,h]  <- A[,,h] - matrix(xi1,dim(A)[1],dim(A)[2])
    B[1,,h] <- B[1,,h] - h * xi0_kappa0

    a[,,h] <- -A[,,h]/h
    b[,,h] <- -B[,,h]/h
  }

  return(list(A=A,B=B,
              a=a,b=b))
}

compute_AB_thk <- function(xi0,xi1,
                           u,
                           H, # maximum maturity,
                           k, # indexation lag of the payoff
                           psi, psi.parameterization){
  # This function is used to price options.
  # It employs recursive formulas to price payoffs of the type:
  #    exp(u'w_{t+h-k}), settled on date t+h.
  # The risk-free short-term rate is r_t = xi0 + xi1'w_t,
  #    where w_t is the state vector, of dimension n x 1.
  # The function returns a list containing A, B, a, and b.
  # Bond prices are such that the date-t payoff is = exp(B_h + A_h'w_t), where
  #    B_h is a scalar and A_h is a (n x 1) vector.
  # u can be of dimension n x q with q > 1 (n is the dimension of the state)
  # psi and psi.parameterization defined the risk-neutral dynamics of w_t.
  # The function returns NaN for horizons h < k.

  n <- length(xi1) # dimension of state vector
  u <- matrix(u,nrow=n) # make sure u has the right dimension

  A.h_1 <- 0
  B.h_1 <- 0

  n <- dim(u)[1]
  q <- dim(u)[2]

  A <- array(NaN, c(n, q, H))
  B <- array(NaN, c(1, q, H))
  for (h in 1:H) {

    X <- matrix(-xi1, n, q)*(h>1) + A.h_1

    psi.u <- psi(X, psi.parameterization)
    A.h <- matrix(psi.u$a, n, q) + u*(h==k)
    B.h <- psi.u$b + B.h_1
    A[, , h] <- A.h
    B[, , h] <- B.h
    A.h_1 <- A.h
    B.h_1 <- B.h
  }

  for(h in 1:H){
    A[,,h]  <- A[,,h] - matrix(xi1, n, q) # one xi1 is missing
    B[1,,h] <- B[1,,h] - h * xi0          # all xi0 are missing
  }

  if(k>1){
    A[,,1:(k-1)] <- NaN
    B[,,1:(k-1)] <- NaN
  }

  return(list(A=A,B=B))
}



# loglik_Gaussian <- function(std_dev,epsilon){
#   # epsilon is of dimension T x n, and epsilon_t ~ N(0,Omega),
#   # where Omega = diag(std_dev^2).
#   # the function returns the log-likelihood associated with those
#   # observations included in matrix epsilon.
#   # epsilon is of dimension T x n
#   # std_dev is of dimension T x n
#   # !!!: works only when components of epsilon are independent
#
#   if(sum(dim(std_dev)==dim(epsilon))!=2){
#     print("Error in loglik_Gaussian: std_dev and epsilon should have the same dimension.")
#     return(NaN)
#   }
#
#   epsilon_standardized <- epsilon * (1/std_dev)
#   logl <- - log(abs(std_dev)) - .5 * epsilon_standardized^2
#
#   return(sum(logl))
# }

log_mvdnorm <- function(X, Sigma){
  # This function computes the log-likelihood associated with a matrix of observations X.
  # The rows of X are supposed to be i.i.d., drawn in N(0,Sigma), where
  # Sigma is a covariance matrix.

  TT <- nrow(X)
  n <- ncol(X)

  Sigma_inv <- solve(Sigma)
  det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus

  term1 <- -0.5 * TT * n * log(2 * pi)
  term2 <- -0.5 * TT * det_Sigma
  term3 <- -0.5 * sum(((X %x% matrix(1,1,n))*(matrix(1,1,n) %x% X)) %*% c(Sigma_inv))

  logL <- term1 + term2 + term3
  return(as.numeric(logL))
}


#
# model <- list(alpha = matrix(c(0,.02),ncol=1),
#               nu = matrix(c(.1,.2),ncol=1),
#               mu = matrix(c(1,1),ncol=1),
#               beta = matrix(c(.9,0,.2,.8),2,2),
#               n_w=2)
# VAR_representation <- compute_expect_variance(psi.VARG,model)
#
# theta1 <- matrix(rnorm(8),4,2)
# theta2 <- matrix(rnorm(8),4,2)
# resEV <- compute_expect_variance_H(VAR_representation,
#                                    theta1=theta1,
#                                    theta2=theta2,H=H)
#
# w0 <- VAR_representation$Ew
#
# N <- 100000
# H <- 5
# all_w <- array(NaN,c(2,N,H+1))
# for(i in 1:N){
#   simres <- simul_VARG(model,nb_periods=5,w0 = w0)
#   all_w[,i,] <- simres
# }
#
# X <- theta1 %*% all_w[,,H+1]
# for(j in 1:(H-1)){
#   X <- X + theta2 %*% all_w[,,H+1-j]
# }
#
# cbind(
#   resEV$all_C0[,,H] + resEV$all_C1[,,H] %*% w0,
#   apply(X,1,mean)
# )
#
# cbind(
#   resEV$all_Gamma0[,,H] + resEV$all_Gamma1[,,H] %*% w0,
#   c(var(t(X)))
# )

price_IR_caps_floors <- function(W, # Values of state vector (T x n)
                                 H, # maximum maturity, in model periods
                                 tau, # maturity of ref. rate, in model periods
                                 freq, # number of model periods per year
                                 all_K, # vector of (annualized) strikes
                                 psi, # Laplace transform of W
                                 parameterization, # see details below
                                 max_x = 100000, # settings for Riemann sum comput.
                                 dx_statio = 10,
                                 min_dx = 1e-05,
                                 nb_x1 = 5000
){
  # W is the matrix of state-vector values, of dimension T x n (T = number of dates)
  # "parameterization" includes:
  #     - "model" (model should include "n_w", the dimension of w_t)
  #     - "xi0", "xi1" (parameterization of short-term rate, not annualized)

  xi0 <- parameterization$xi0
  xi1 <- parameterization$xi1

  alpha_tau <- tau/freq # year fraction corresponding to tau

  varphi <- function(x,parameterization,H){
    # This function is an argument of truncated.payoff

    u <- parameterization$u
    v <- parameterization$v # determine thresholds (dimension nw x 1)
    i.v.x <- matrix(v,ncol=1) %*% matrix(c(1i*x),nrow=1)

    tau <- parameterization$tau

    nw <- dim(i.v.x)[1]
    q  <- dim(i.v.x)[2]
    U <- matrix(u, nw, q) + i.v.x

    res <- compute_AB_thk(xi0 = parameterization$xi0,
                          xi1 = parameterization$xi1,
                          u = U,
                          H = H, # maximum maturity, in model periods
                          k = tau, # indexation lag of the payoff
                          psi = psi,
                          psi.parameterization = parameterization$model)
    return(res)
  }

  # Compute affine representation of the reference rate:
  res_tau <- compute_AB_classical(xi0 = parameterization$xi0,
                                  xi1 = parameterization$xi1,
                                  H = tau,
                                  psi = psi,
                                  psi.parameterization = parameterization$model)
  # Specification of discount factor:
  A_tau <- res_tau$A[,,tau]
  B_tau <- res_tau$B[,,tau]
  a_tau <- res_tau$a[,,tau]
  b_tau <- res_tau$b[,,tau]

  thresholds <- - log(1 + alpha_tau*all_K) - B_tau
  b.matrix <- t(matrix(thresholds,length(thresholds),H))

  TT <- dim(W)[1] # number of dates.
  vec1TT <- matrix(1,TT,1)

  parameterization$tau <- tau
  parameterization$u   <- NaN

  # Use of "truncated.payoff" for caplets:
  parameterization$v <- matrix(+ A_tau,ncol=1)
  parameterization$u <- matrix(- A_tau,ncol=1)
  res_truncated1 <- truncated.payoff(W,
                                     b.matrix = b.matrix,
                                     varphi = varphi,
                                     parameterization = parameterization,
                                     max_x = max_x,
                                     dx_statio = dx_statio,
                                     min_dx = min_dx,
                                     nb_x1 = nb_x1)
  parameterization$u <- matrix(0*A_tau, ncol = 1)
  res_truncated2 <- truncated.payoff(W,
                                     b.matrix = b.matrix,
                                     varphi = varphi,
                                     parameterization = parameterization,
                                     max_x = max_x,
                                     dx_statio = dx_statio,
                                     min_dx = min_dx,
                                     nb_x1 = nb_x1)

  # Use of "truncated.payoff" for floorlets:
  b.matrix <- - t(matrix(thresholds,length(thresholds),H))
  parameterization$v <- matrix(- A_tau,ncol=1)
  parameterization$u <- matrix(0*A_tau,ncol=1)
  res_truncated3 <- truncated.payoff(W,
                                     b.matrix = b.matrix,
                                     varphi = varphi,
                                     parameterization = parameterization,
                                     max_x = max_x,
                                     dx_statio = dx_statio,
                                     min_dx = min_dx,
                                     nb_x1 = nb_x1)
  parameterization$u <- matrix(- A_tau, ncol = 1)
  res_truncated4 <- truncated.payoff(W,
                                     b.matrix = b.matrix,
                                     varphi = varphi,
                                     parameterization = parameterization,
                                     max_x = max_x,
                                     dx_statio = dx_statio,
                                     min_dx = min_dx,
                                     nb_x1 = nb_x1)

  matrix_K <- t(matrix(all_K,length(all_K),TT))
  array_K  <- array(matrix_K,c(TT,length(all_K),H))

  Caplet <- exp(-B_tau) * res_truncated1 -
    (1 + alpha_tau * array_K) * res_truncated2

  Floorlet <- (1 + alpha_tau * array_K) * res_truncated3 -
    exp(-B_tau) * res_truncated4


  # Compute maturities available given H and tau:
  maturities.in.model.period <- seq(tau,H,by=tau)

  Caps <- array(NaN,c(TT,
                      length(all_K),
                      length(maturities.in.model.period)))
  dimnames(Caps) <- list(
    time     = paste0("t", 1:TT),
    strike   = paste0(100*all_K,"%"),
    maturity = paste0("matur_",maturities.in.model.period/freq,"_years")
  )
  Floors <- Caps

  Cap   <- 0
  Floor <- 0
  count <- 0
  for(maturity in maturities.in.model.period){
    count <- count + 1
    Cap   <- Cap   + Caplet[,,maturity]
    Floor <- Floor + Floorlet[,,maturity]
    Caps[,,count]   <- Cap
    Floors[,,count] <- Floor
  }

  return(
    list(A_tau = A_tau, B_tau = B_tau,
         a_tau = a_tau, b_tau = b_tau,
         Caplet = Caplet, Floorlet = Floorlet,
         Caps = Caps,     Floors = Floors,
         maturities.in.model.period = maturities.in.model.period)
  )
}



