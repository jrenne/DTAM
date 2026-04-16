#' Fast Gaussian multi-horizon aggregation
#'
#' Computes the mean, autoregressive loading, and covariance associated with the
#' cumulated Gaussian VAR state `x_{t+1} + ... + x_{t+h}` for a collection of
#' horizons.
#'
#' @param model List containing `mu`, `Phi`, and `Sigma` for the Gaussian
#'   VAR(1).
#' @param Maturities_decompo Integer decomposition matrix describing the target
#'   maturities.
#'
#' @return A list with arrays `MU`, `PHI`, and `SIGMA`, plus a vector `H` of
#'   maturities.
#'
#' @examples
#' model <- list(
#'   mu = matrix(c(0, 0.1), 2, 1),
#'   Phi = matrix(c(0.8, 0.1,
#'                  0.0, 0.7), 2, 2, byrow = TRUE),
#'   Sigma = diag(c(0.2, 0.1))
#' )
#' decomp <- matrix(c(2, 0,
#'                    1, 2), 2, 2, byrow = TRUE)
#' res <- compute_Gaussian_fast(model, decomp)
#' res$H
#'
#' @export
compute_Gaussian_fast <- function(model,Maturities_decompo){
  # This function computes mu_h, Phi_h, and Sigma_h that are s.t.
  # x_{t+h}+...+x_{t+1}|N(mu_h + Phi_h·x_t, Sigma_h),
  # where x_t follows the Gaussian VAR(1):
  # x_{t+1} = mu + Phi·x_t + Sigma^{1/2}·epsilon_{t+1},
  #          where epsilon_{t+1}~N(O,Id).
  # Maturities_decompo determines (and decomposes)
  # these maturities; it is of dimension k x k. The k-th row determines
  # how the k-th horizon will be decomposed into the previous ones. Denote
  # the (i,j) entry of Maturities_decompo with m(i,j). Denote the i-th horizon
  # with h(i). We have:
  # h(i) = m(i,1)*1 + m(i,2)*h(1) + ... + m(i,i-1)*h(i-1).
  # For instance:
  #   Maturities_decompo[1,] = [2, 0, 0] => h_1 = 2
  #   Maturities_decompo[2,] = [1, 2, 0] => h_2 = 1*1 + 2*2 = 5
  #   Maturities_decompo[3,] = [2, 0, 2] => h_3 = 2*1 + 0*2 + 2*5 = 12.
  # This decomposition speeds up the calculation by decreasing the number of
  # matrix multiplications.

  mu  <- model$mu
  Phi <- model$Phi
  Sigma <- model$Sigma # covariance matrix

  n <- dim(Phi)[1]

  k <- dim(Maturities_decompo)[1]

  Id_n  <- diag(n)
  Id_n2 <- diag(n*n)
  Id_Phi_1 <- solve(Id_n - Phi)

  PhiPhi <- Phi %x% Phi
  Id_PhiPhi_1 <- solve(Id_n2 - PhiPhi)

  Sigma_star <- Id_Phi_1 %*% Sigma %*% t(Id_Phi_1)

  MU    <- array(NaN,c(n,1,k))
  PHI   <- array(NaN,c(n,n,k))
  SIGMA <- array(NaN,c(n,n,k))
  H <- NULL

  all_Phi_h <- array(NaN,c(n,n,k))
  all_PhiPhi_h <- array(NaN,c(n*n,n*n,k))

  for(i in 1:k){
    Phi_exp_h    <- Id_n
    PhiPhi_exp_h <- Id_n2
    h_i <- 0
    for(j in 1:i){
      if(j==1){
        Aux <- Phi
        AAux <- PhiPhi
        aux <- 1
      }else{
        Aux <- all_Phi_h[,,j-1]
        AAux <- all_PhiPhi_h[,,j-1]
        aux <- H[j-1]
      }
      Phi_exp_h    <- Phi_exp_h    %*% (Aux %^% Maturities_decompo[i,j])
      PhiPhi_exp_h <- PhiPhi_exp_h %*% (AAux %^% Maturities_decompo[i,j])
      h_i  <- h_i + aux*Maturities_decompo[i,j]
    }
    all_Phi_h[,,i]    <- Phi_exp_h
    all_PhiPhi_h[,,i] <- PhiPhi_exp_h

    K_h   <- Id_Phi_1 %*% (Id_n - Phi_exp_h)
    mu_h  <- Id_Phi_1 %*% (h_i*Id_n - Phi %*% K_h) %*% mu
    Phi_h <- Phi %*% K_h
    aux <- Phi %*% K_h %*% Sigma_star

    Aux <- PhiPhi %*% Id_PhiPhi_1 %*% (Id_n2 - PhiPhi_exp_h) %*% c(Sigma_star)
    vecSigma_h <- h_i*c(Sigma_star) - c(aux) - c(t(aux)) + Aux

    # Store results:
    MU[,,i]  <- mu_h
    PHI[,,i] <- Phi_h
    SIGMA[,,i] <- matrix(vecSigma_h,n,n)
    H <- c(H,h_i)
  }
  return(list(MU=MU,PHI=PHI,SIGMA=SIGMA,H=H))
}

# Maturities_decompo <- diag(4)
# Maturities_decompo[2,1:2] <- c(0,4)
# Maturities_decompo[3,1:3] <- c(0,1,3)
# Maturities_decompo[4,1:4] <- c(1,1,1,4)
#
# mu <- matrix(c(0,1),2,1)
# Phi <- .9*diag(2)
# Phi[2,1] <- .4
# Sigma <- .1*diag(2)
# Sigma[1,2] <- .005
# Sigma[2,1] <- .005
#
# model <- list(mu=mu,Phi=Phi,Sigma=Sigma)
#
# RES <- compute_Gaussian_fast(model,Maturities_decompo)
# maxH <- RES$H[length(RES$H)]
#
# u <- matrix(1,2,1)
# res <- reverse.MHLT(psi.GaussianVAR,u1 = u,H = maxH,psi.parameterization = model)
#
# m <- dim(Maturities_decompo)[1]
# mu_h <- matrix(RES$MU[,,m],ncol=1)
# Phi_h <- RES$PHI[,,m]
# Sigma_h <- RES$SIGMA[,,m]
#
# t(u) %*% mu_h + (t(u) %*% Sigma_h %*% u)/2
# t(u) %*% Phi_h

simul.GVAR_OLD <- function(model,nb.sim,x0=NaN){
  n <- dim(model$Phi)[1]
  if(is.na(x0[1])){
    x0 <- solve(diag(n) - model$Phi) %*% model$mu
  }

  if(is.null(model$Sigma12)){
    Sigma12 <- t(chol(model$Sigma))
  }else{
    Sigma12 <- model$Sigma12
  }

  X <- c(x0)
  x <- x0
  for(t in 2:nb.sim){
    x <- model$mu + model$Phi %*% x + Sigma12 %*% rnorm(dim(Sigma12)[2])
    X <- rbind(X,c(x))
  }
  return(X)
}

#' Simulate a Gaussian VAR process
#'
#' Simulates trajectories from a Gaussian VAR(1) model of the form
#' \deqn{
#' x_t = \mu + \Phi x_{t-1} + \Sigma^{1/2}\varepsilon_t,
#' }{
#' x_t = mu + Phi x_{t-1} + Sigma^{1/2} eps_t,
#' }
#' where \eqn{\varepsilon_t \sim \mathcal{N}(0, I)}{eps_t ~ N(0, I)}.
#'
#' @param model A list describing the VAR. It must contain `mu` and `Phi`, and
#'   either `Sigma12` (the square-root innovation matrix) or `Sigma`
#'   (the innovation covariance matrix).
#' @param nb.sim Number of simulated dates.
#' @param x0 Initial state vector. If left at the default, the function starts
#'   from the unconditional mean implied by `mu` and `Phi`.
#' @param nb.replic Number of independent simulated paths. When `nb.replic = 1`,
#'   the output is a matrix. For `nb.replic > 1`, the output is a three-way
#'   array.
#'
#' @return If `nb.replic = 1`, a matrix of dimension `nb.sim x n`, where `n` is
#'   the dimension of the state vector. If `nb.replic > 1`, an array of
#'   dimension `nb.sim x n x nb.replic`.
#'
#' @details
#' If `model$Sigma12` is not provided, the function computes it as
#' `t(chol(model$Sigma))`.
#'
#' This function is used repeatedly in the Bookdown companion project to
#' generate sample paths for affine Gaussian models before pricing or
#' estimation steps.
#'
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#'
#' @examples
#' # Example adapted from the Bookdown companion project.
#' set.seed(123)
#'
#' model <- list(
#'   mu = matrix(c(0, 0), ncol = 1),
#'   Phi = matrix(c(0.95, 0.05,
#'                  0.00, 0.90), nrow = 2, byrow = TRUE),
#'   Sigma12 = diag(c(0.1, 0.2))
#' )
#'
#' sim <- simul.GVAR(model, nb.sim = 50)
#' dim(sim)
#'
#' sim_multi <- simul.GVAR(model, nb.sim = 20, nb.replic = 3)
#' dim(sim_multi)
#'
#' @export
simul.GVAR <- function(model,nb.sim,x0=NaN,nb.replic=1){
  # Simulate a Gaussian VAR model.
  # If nb.replic > 1, simulate several processes of length nb.sim in parallel;
  # the output then is an array (three dimensions).

  n <- dim(model$Phi)[1]
  if(is.na(x0[1])){
    x0 <- solve(diag(n) - model$Phi) %*% model$mu
  }

  if(is.null(model$Sigma12)){
    Sigma12 <- t(chol(model$Sigma))
  }else{
    Sigma12 <- model$Sigma12
  }

  if(nb.replic==1){
    X <- c(x0)
    x <- x0
    for(t in 2:nb.sim){
      x <- model$mu + model$Phi %*% x + Sigma12 %*% rnorm(dim(Sigma12)[2])
      X <- rbind(X,c(x))
    }
  }else{
    # In that case, we simulate nb.replic paths in parallel
    X <- array(NaN,c(nb.sim,n,nb.replic))
    X[1,,] <- x0
    x <- matrix(x0,n,nb.replic)
    for(t in 2:nb.sim){
      x <- matrix(model$mu,n,nb.replic) + model$Phi %*% x + Sigma12 %*% matrix(rnorm(dim(Sigma12)[2]*nb.replic),dim(Sigma12)[2],nb.replic)
      X[t,,] <- x
    }
  }
  return(X)
}
#' Plot Gaussian iso-density contours
#'
#' Draws contour lines associated with a bivariate Gaussian distribution.
#'
#' @param Ew Mean vector of length two.
#' @param Vw `2 x 2` covariance matrix.
#' @param n Grid size used for contour evaluation.
#' @param min_w,max_w Optional plotting bounds. If omitted, the routine uses a
#'   multiple of the marginal standard deviations.
#' @param prob_levels Probability-content levels used to determine contour
#'   radii.
#' @param n_std Number of standard deviations used when bounds are not supplied.
#' @param maintitle Plot title.
#'
#' @return Invisibly returns a list containing the plotting grid and contour
#'   levels.
#'
#' @examples
#' plot_isodensity_gaussian(Ew = c(0, 0),
#'                          Vw = matrix(c(1, 0.5, 0.5, 2), 2, 2),
#'                          maintitle = "")
#'
#' @export
