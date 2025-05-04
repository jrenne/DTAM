
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
  # psi admits 2 arguments: model and u; model is a list of parameters; u is the
  # vector of points where one wants to evaluate the pdf.

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


compute_H_bar_TopDown <- function(model,gamma,H=10,W=NaN){
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

  varphi <- reverse.MHLT(psi.w.TopDown,
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


compute_CDS_TopDown <- function(model,H=10,W,RR=.5){

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
    res <- compute_H_bar_TopDown(model,ej_tilde,H=10,W)
    numerator <- 0
    denominat <- 0
    for(h in 1:H){
      numerator <- numerator + res$H_upper[,h] - res$H_lower[,h]
      denominat <- denominat + I_j * res$B[,h]       - res$H_lower[,h]
      CDS[,j,h] <- (1-RR)*(numerator/denominat)
    }
  }

  return(CDS)
}

varphi4G_TopDown <- function(u,parameterization){

  model   <- parameterization$model
  H       <- parameterization$H # maturity
  gamma   <- parameterization$gamma
  indic_Q     <- parameterization$indic_Q # if FALSE, use physical LT, otherwise Q
  indic_upper <- parameterization$indic_upper # if TRUE, compute G_upper, otherwise G_lower

  nw <- dim(u)[1]
  q  <- dim(u)[2]

  u2 <- matrix(gamma - model$xi1, nw, q)
  u1 <- u
  if(indic_upper){
    u1    <- u1 + matrix(gamma, nw, q)
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

  return(list(A = res_reverse$A,
              B = matrix(res_reverse$B,nrow=1)))
}

truncated.payoff <- function(model,
                             W, # values of w_t
                             gamma,v,vector.of.b,
                             H,
                             varphi,
                             parameterization,
                             max_x = 500,
                             dx_statio = .1,
                             min_dx = 1e-05,
                             nb_x1 = 1000){
  # This function calculates the value of payoff(w)·1_{v'w < b}, for different
  #      values of b (collected in "vector.of.b").
  #       "W" is of dimension T x n.w; it contains values of the state vector w
  #          at which we want to evaluate the prices.
  #       "varphi" is a function that takes two arguments:
  #          (i) u and (ii) parameterization;
  #       it returns matrices A(u) and B(u) that are such that the "price" of
  #          {payoff x exp(u)} is given by exp(A(u)'w + B(u)),
  #          where A(u) is of dimension n.u x n.w,
  #          and B(u) is of dimension n.u x 1
  # The truncation is of the form v'w < b,
  #          "vector.of.b" collects the b boundaries of interest.

  n_w <- model$n_w

  # x <- matrix(exp(seq(log(min_x),log(max_x),length.out=nb_x)),ncol=1)
  # x <- matrix(seq(min_x,max_x,length.out=nb_x),ncol=1)
  #
  # x1 <- matrix(exp(seq(log(min_x),log(max_x_log),length.out=nb_x/2)),ncol=1)
  # x2 <- matrix(seq(max_x_log,max_x,length.out=nb_x/2),ncol=1)
  # x <- rbind(x1,x2)

  dx1 <- exp(seq(log(min_dx),log(dx_statio),length.out=nb_x1))
  max_x1 <- tail(cumsum(dx1),1)
  nb_x2 <- (max_x - max_x1)/dx_statio
  dx <- c(dx1,rep(dx_statio,nb_x2))
  x <- matrix(cumsum(dx),ncol=1)

  nb_x <- length(x)

  i.v.x <- matrix(v,ncol=1) %*% matrix(c(1i*x),nrow=1)

  # Add necessary fields in parameterization:
  parameterization$model <- model
  parameterization$gamma <- gamma
  parameterization$H     <- H

  res_varphi  <- varphi(i.v.x,parameterization)
  res_varphi0 <- varphi(matrix(0,n_w,1),parameterization)

  xi0 <- model$xi0
  xi1 <- model$xi1

  B <- matrix(res_varphi$B,1,nb_x*H) -
    matrix(1:H,nrow=1) %x% matrix(1,1,nb_x) * xi0
  A <- matrix(res_varphi$A,n_w,nb_x*H) - matrix(xi1,n_w,nb_x*H)

  B0 <- matrix(res_varphi0$B,nrow=1) - (1:H) * xi0
  A0 <- matrix(res_varphi0$A,n_w,H)  - matrix(xi1,n_w,H)

  if(dim(W)[2]!=n_w){# to make sure W is of dimension T x n.w
    W <- t(W)
  }
  T <- dim(W)[1]

  psi_eval  <- exp(t(A)  %*% t(W) + matrix(B, nb_x*H,T))
  psi_eval0 <- exp(t(A0) %*% t(W) + matrix(B0,H,T))

  # modify format to have them nb_x x (T*H):
  Psi_eval  <- matrix(NaN,nb_x,T*H)
  for(h in 1:H){
    Psi_eval[,(T*(h-1)+1):(T*h)]  <- psi_eval[(nb_x*(h-1)+1):(nb_x*h),]
  }

  dx <- matrix(x - c(0,x[1:(nb_x-1)]),nb_x,T*H)
  x_mid <- matrix(x,nb_x,T*H)

  all_prices <- matrix(0,T*H,length(vector.of.b))

  count_b <- 0
  for(b in vector.of.b){
    count_b <- count_b + 1
    aux <- Im(Psi_eval * matrix(exp(-1i*x_mid*b),nb_x,T*H))
    fx <- aux/x_mid * dx
    f  <- 1/2*c(t(psi_eval0)) - 1/pi * apply(fx,2,sum)

    all_prices[,count_b] <- f
  }

  All_prices <- array(NaN,c(T,length(vector.of.b),H))
  for(h in 1:H){
    All_prices[,,h] <- all_prices[(T*(h-1)+1):(T*h),]
  }

  return(All_prices)
}

compute_G <- function(model,W,
                      gamma,v,
                      vector.of.b, # attachment/detachment points
                      H, # maturity
                      indic_Q = TRUE, # computed under Q.
                      indic_upper = TRUE, # otherwise, computed G_lower
                      max_x = 2000,
                      dx_statio = 2,
                      min_dx = 10^(-6),
                      nb_x1=1000){
  # This computes functions G_lower and G_upper

  parameterization <- list(indic_Q = indic_Q,
                           indic_upper = indic_upper)

  P <- truncated.payoff(model,W,
                        gamma,v,vector.of.b,
                        H,
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
  # This function computes the price of tranche products.
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
                            vector.of.b=vector.of.a.bar, # attachment/detachment points
                            H=H,
                            indic_Q = indic_Q,
                            indic_upper = TRUE,
                            max_x = max_x,dx_statio = dx_statio,
                            min_dx = min_dx,nb_x1 = nb_x1)
  G_lower_0_ej <- compute_G(model,
                            W=W,
                            gamma = 0*ej_tilde,
                            v = ej_tilde,
                            vector.of.b=vector.of.a.bar, # attachment/detachment points
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
                              vector.of.b=vector.of.a.bar, # attachment/detachment points
                              H=H,
                              indic_Q = indic_Q,
                              indic_upper = TRUE,
                              max_x = max_x,dx_statio = dx_statio,
                              min_dx = min_dx,nb_x1 = nb_x1)
  G_lower_uej_ej <- compute_G(model,
                              W=W,
                              gamma = u*ej_tilde,
                              v = ej_tilde,
                              vector.of.b=vector.of.a.bar, # attachment/detachment points
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


