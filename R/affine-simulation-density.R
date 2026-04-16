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

#' Simulate scalar affine count and intensity processes
#'
#' `simul.ARG()` simulates a scalar ARG/ARG0 process. `simul.compound.poisson()`
#' simulates the compound-Poisson process used in the affine-process chapter.
#'
#' @param nb.sim Number of simulated dates.
#' @param mu,nu,rho,alpha Scalar ARG parameters.
#' @param w0 Optional initial value.
#' @param nb.replic Number of independent ARG paths simulated in parallel.
#' @param Gamma,Pi,lambda Parameters of the compound-Poisson process.
#'
#' @return `simul.ARG()` returns a `nb.sim x nb.replic` matrix. When
#'   `nb.replic = 1`, this is a single-column matrix. `simul.compound.poisson()`
#'   returns a numeric vector of length `nb.sim`.
#'
#' @examples
#' set.seed(123)
#' sim_arg <- simul.ARG(nb.sim = 20, mu = 0.05, nu = 1, rho = 0.95, alpha = 0)
#' dim(sim_arg)
#'
#' sim_cp <- simul.compound.poisson(nb.sim = 20, Gamma = 1, Pi = 0.95,
#'                                  lambda = 0.1)
#' length(sim_cp)
#'
#' @name simul_ARG
#' @export
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


#' @rdname simul_ARG
#' @export
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

#' Recover a density from a Laplace transform
#'
#' Uses Fourier inversion to approximate a univariate density on a user-supplied
#' grid of values.
#'
#' @param model List of model parameters passed to `psi`.
#' @param values.of.variable Grid at which the density is approximated.
#' @param psi Function returning the Laplace transform at the requested points.
#' @param max_x,min_x,nb_x Numerical integration settings for the inversion.
#'
#' @return Numeric vector of approximate density values on
#'   `values.of.variable`.
#'
#' @examples
#' psi_exp <- function(u, model) {
#'   -log(1 - model$mu * u)
#' }
#' grid <- seq(0.01, 2, length.out = 50)
#' dens <- make.pdf(list(mu = 1), values.of.variable = grid, psi = psi_exp,
#'                  max_x = 100, nb_x = 500)
#' length(dens)
#'
#' @export
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

#' Compute conditional and unconditional moments from an affine transform
#'
#' Recovers the affine conditional mean and variance representation of a state
#' vector from its conditional log-Laplace transform, and deduces the
#' corresponding unconditional moments.
#'
#' @param psi Conditional log-Laplace transform of the state vector. It must
#'   accept arguments `(u, model)` and return a list with components `a` and
#'   `b`.
#' @param model A list containing the parameterization of the state process. It
#'   must include `n_w`, the dimension of the state vector, together with the
#'   objects required by `psi`.
#' @param du Finite-difference step used to approximate first and second
#'   derivatives of the log-Laplace transform around zero.
#'
#' @return A list with:
#'   `mu`, `Phi`, `Theta0`, `Theta1` describing the affine conditional moments,
#'   and `Ew`, `Vw` giving the unconditional mean and variance.
#'
#' @details
#' The function numerically differentiates the conditional log-Laplace
#' transform:
#'
#' \deqn{
#' \mathbb{E}_t[\exp(u^\prime w_{t+1})] =
#' \exp\left(a(u)^\prime w_t + b(u)\right),
#' }{
#' E_t[exp(u' w_{t+1})] = exp(a(u)' w_t + b(u)),
#' }
#'
#' around `u = 0` to recover:
#'
#' \deqn{
#' \mathbb{E}_t(w_{t+1}) = \mu + \Phi w_t,
#' }{
#' E_t(w_{t+1}) = mu + Phi w_t,
#' }
#'
#' and
#'
#' \deqn{
#' \mathrm{vec}\{\mathrm{Var}_t(w_{t+1})\} = \Theta_0 + \Theta_1 w_t.
#' }{
#' vec(Var_t(w_{t+1})) = Theta0 + Theta1 w_t.
#' }
#'
#' These objects are used throughout the package for recursive preferences,
#' shadow-rate approximations, and filtering applications.
#'
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#'
#' @examples
#' model <- list(
#'   mu = matrix(c(0.01, 0.02), ncol = 1),
#'   Phi = matrix(c(0.95, 0.00,
#'                  0.05, 0.90), nrow = 2, byrow = TRUE),
#'   Sigma = diag(c(0.02, 0.01)),
#'   n_w = 2
#' )
#'
#' res <- compute_expect_variance(psi.GaussianVAR, model)
#'
#' res$mu
#' res$Phi
#' res$Ew
#'
#' @export
