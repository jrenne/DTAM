#' Rating-migration bond pricing and simulation
#'
#' `compute_bond_price_Ratings()` prices risky and risk-free bonds in a model
#' with stochastic rating migration. `simul.rating.migration()` simulates a
#' rating path conditional on an exogenous state trajectory.
#'
#' @param model List of model parameters.
#' @param Y Matrix of factor realizations.
#' @param H Maximum maturity.
#' @param psi One-step Laplace transform of the factor process.
#' @param tau_ini Initial rating for the simulated migration path.
#'
#' @return `compute_bond_price_Ratings()` returns a list containing risk-free and
#'   risky bond prices, yields, and spreads. `simul.rating.migration()` returns
#'   a vector of simulated ratings.
#'
#' @examples
#' model <- list(
#'   kappa0 = matrix(c(0.01, 0), 1, 2),
#'   kappa1 = matrix(c(0.1, 0), 1, 2),
#'   xi0 = 0.01,
#'   xi1 = matrix(0.1, 1, 1),
#'   V = diag(2),
#'   psi.parameterization = list(mu = matrix(0, 1, 1),
#'                               Phi = matrix(0.8, 1, 1),
#'                               Sigma12 = matrix(0.1, 1, 1))
#' )
#' Y <- matrix(c(0, 0.1, -0.1), 3, 1)
#' res <- compute_bond_price_Ratings(model, Y, H = 2, psi = psi.GaussianVAR)
#' dim(res$RF_bond_price)
#' simul.rating.migration(model, Y, tau_ini = 1)
#'
#' @name Ratings_tools
#' @export
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


#' @rdname Ratings_tools
#' @export
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


#' Compute affine bond-pricing coefficients
#'
#' Computes the affine coefficients used to price risk-free bonds and, when
#' default-intensity parameters are supplied, defaultable bonds.
#'
#' Bond prices take the form
#' \deqn{
#' B_{t,h} = \exp\left(B_h + A_h^\prime w_t\right),
#' }{
#' B_{t,h} = exp(B_h + A_h' w_t),
#' }
#' and yields are
#' \deqn{
#' y_{t,h} = b_h + a_h^\prime w_t.
#' }{
#' y_{t,h} = b_h + a_h' w_t.
#' }
#'
#' @param xi0 Scalar intercept in the short-rate equation.
#' @param xi1 Loading vector in the short-rate equation.
#' @param kappa0 Optional intercept of the default intensity. If left as `NaN`,
#'   only risk-free bonds are priced.
#' @param kappa1 Optional loading matrix for default intensities. Different
#'   columns can correspond to different entities.
#' @param H Maximum maturity, expressed in model periods.
#' @param psi Conditional Laplace transform of the state vector.
#' @param psi.parameterization List of parameters passed to `psi`.
#'
#' @return A list with four arrays: `A`, `B`, `a`, and `b`. The first two give
#'   the affine price coefficients, and the latter two the corresponding yield
#'   coefficients.
#'
#' @details
#' If `kappa0` and `kappa1` are provided, the function computes both risk-free
#' and defaultable bond coefficients. Otherwise it returns only the risk-free
#' coefficients.
#'
#' Internally, the recursion relies on `reverse.MHLT()`.
#'
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#'
#' @examples
#' # Example adapted from the Bookdown companion project.
#' model <- list(
#'   mu = matrix(c(0.01, 0.02), ncol = 1),
#'   Phi = matrix(c(0.95, 0.00,
#'                  0.05, 0.90), nrow = 2, byrow = TRUE),
#'   Sigma = diag(c(0.02, 0.01))
#' )
#'
#' xi0 <- 0
#' xi1 <- matrix(c(1, 0), ncol = 1)
#'
#' resAB <- compute_AB_classical(
#'   xi0 = xi0,
#'   xi1 = xi1,
#'   H = 12,
#'   psi = psi.GaussianVAR,
#'   psi.parameterization = model
#' )
#'
#' dim(resAB$A)
#' dim(resAB$a)
#'
#' @export
