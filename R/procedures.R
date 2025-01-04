
reverse.MHLT <- function(psi,u1,u2=NaN,H,psi.parameterization){
  # compute the multi-horizon Laplace transform in the reverse-order case
  # That is, we consider:
  # E_t(exp(u_2'w_{t+1}+...+u_2'w_{t+h-1}+u_1'w_{t+h}))
  # Inputs: psi is the (one-period) conditional Laplace transform of process w_t
  # The function can evaluate the MHLT in parallel for k different vectors
  #    u1 and u2 (i.e., u1 and u2, as inputs can be of dimension n x k with k > 1)
  if(is.na(u2)){# in that case, u2 <- u1
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
  # If w_t is n-dimensional, u is of dimension n x k
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
