#' Price interest-rate caps and floors
#'
#' Prices caplets/floorlets on an affine reference rate and accumulates them
#' into caps and floors.
#'
#' @param W Matrix of state-vector values.
#' @param H Maximum maturity in model periods.
#' @param tau Maturity of the reference rate, in model periods.
#' @param freq Number of model periods per year.
#' @param all_K Vector of annualized strikes.
#' @param psi One-step Laplace transform of the state vector.
#' @param parameterization List containing the model and the short-rate
#'   specification.
#' @param max_x,dx_statio,min_dx,nb_x1 Numerical integration settings.
#'
#' @return A list containing `Caplet`, `Floorlet`, `Caps`, `Floors`, and the
#'   affine coefficients of the reference rate.
#'
#' @examples
#' model <- list(
#'   n_w = 1,
#'   mu = matrix(0, 1, 1),
#'   Phi = matrix(0.9, 1, 1),
#'   Sigma12 = matrix(0.1, 1, 1)
#' )
#' W <- matrix(c(0, 0.1), 2, 1)
#' res <- price_IR_caps_floors(
#'   W = W, H = 4, tau = 1, freq = 4, all_K = c(0.01, 0.02),
#'   psi = psi.GaussianVAR,
#'   parameterization = list(model = model, xi0 = 0, xi1 = 0)
#' )
#' dim(res$Caps)
#'
#' @export
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

  nb_strikes <- length(all_K)

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

  matrix_K <- t(matrix(all_K,nb_strikes,TT))
  array_K  <- array(matrix_K,c(TT,nb_strikes,H))

  Caplet <- exp(-B_tau) * res_truncated1 -
    (1 + alpha_tau * array_K) * res_truncated2

  Floorlet <- (1 + alpha_tau * array_K) * res_truncated3 -
    exp(-B_tau) * res_truncated4


  # Compute maturities available given H and tau:
  maturities.in.model.period <- seq(tau,H,by=tau)

  Caps <- array(NaN,c(TT,
                      nb_strikes,
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


#' Price stock options from affine state dynamics
#'
#' Prices European stock calls and puts when log returns are affine in the state
#' vector and the payoff is computed through truncated Fourier inversion.
#'
#' @param W Matrix of state-vector values.
#' @param S Vector of ex-dividend stock prices.
#' @param H Maximum maturity in model periods.
#' @param a,b Stock-return loadings.
#' @param K_over_S Vector of strikes expressed as fractions of the spot price.
#' @param psi One-step Laplace transform of the state vector.
#' @param parameterization List containing at least the pricing model and the
#'   short-rate specification.
#' @param max_x,dx_statio,min_dx,nb_x1 Numerical integration settings.
#'
#' @return A list with arrays `Calls` and `Puts`.
#'
#' @examples
#' model <- list(n_w = 1, mu = matrix(0, 1, 1), Phi = matrix(0.9, 1, 1),
#'               Sigma12 = matrix(0.1, 1, 1))
#' W <- matrix(c(0, 0.1), 2, 1)
#' S <- c(100, 102)
#' res <- price_Stock_calls_puts(
#'   W = W, S = S, H = 2, a = 0.1, b = 0,
#'   K_over_S = c(0.95, 1.05),
#'   psi = psi.GaussianVAR,
#'   parameterization = list(model = model, xi0 = 0, xi1 = 0)
#' )
#' dim(res$Calls)
#'
#' @export
price_Stock_calls_puts <- function(W, # Values of state vector (T x n)
                                   S, # ex-dividend stock price (T x 1)
                                   H, # maximum maturity, in model periods
                                   a, b, # specif. of ex-dividend stock returns
                                   K_over_S, # vector of strikes
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
  # Regarding the specification of stock returns:
  #     - "a" is a vector, "b" is a scalar.
  #     - the returns are expressed at the model frequency
  # K_over_S contain strikes, expressed as fractions of the ex-dividend stock price

  xi0 <- parameterization$xi0
  xi1 <- matrix(parameterization$xi1,ncol=1)

  a <- matrix(a,ncol=1)

  nw <- length(xi1)

  varphi <- function(x,parameterization,H){
    # This function is an argument of truncated.payoff

    u <- parameterization$u
    v <- parameterization$v # determine thresholds (dimension nw x 1)
    i.v.x <- matrix(v,ncol=1) %*% matrix(c(1i*x),nrow=1)

    tau <- parameterization$tau

    nw <- dim(i.v.x)[1]
    q  <- dim(i.v.x)[2]
    U <- matrix(u, nw, q) + i.v.x

    res <- reverse.MHLT(psi = psi,
                        u1 = U,
                        u2 = U,
                        H  = H,
                        psi.parameterization = parameterization$model)
    return(res)
  }

  # Computation of thresholds:
  nb_strikes <- length(K_over_S)
  b.matrix <- t(matrix(log(K_over_S),nb_strikes,H)) -
    (b*matrix(1:H,H,nb_strikes))

  TT <- dim(W)[1] # number of dates.
  vec1TT <- matrix(1,TT,1)

  # Use of "truncated.payoff" for payoff involving S:
  parameterization$v <- a
  parameterization$u <- a
  G_0a0a <- truncated.payoff(W,
                             b.matrix = b.matrix,
                             varphi = varphi,
                             parameterization = parameterization,
                             max_x = max_x,
                             dx_statio = dx_statio,
                             min_dx = min_dx,
                             nb_x1 = nb_x1)
  # Use of "truncated.payoff" for payoff involving K:
  parameterization$v <- a
  parameterization$u <- 0*a
  G_000a <- truncated.payoff(W,
                             b.matrix = b.matrix,
                             varphi = varphi,
                             parameterization = parameterization,
                             max_x = max_x,
                             dx_statio = dx_statio,
                             min_dx = min_dx,
                             nb_x1 = nb_x1)

  # Computation of swap:
  res_aa <- reverse.MHLT(psi = psi,
                         u1 = a,
                         u2 = a,
                         H  = H,
                         psi.parameterization = parameterization$model)
  A_aa <- matrix(res_aa$A,nw,H)
  B_aa <- matrix(res_aa$B,H,1)

  res_00 <- reverse.MHLT(psi = psi,
                         u1 = 0*a,
                         u2 = 0*a,
                         H  = H,
                         psi.parameterization = parameterization$model)
  A_00 <- matrix(res_00$A,nw,H)
  B_00 <- matrix(res_00$B,H,1)

  varphi_th_aa <- exp(W %*% A_aa + vec1TT %*% t(B_aa + b * matrix(1:H,H,1)))
  varphi_th_aa <- varphi_th_aa %x% matrix(1,1,nb_strikes)
  varphi_th_00 <- exp(W %*% A_00 + vec1TT %*% t(B_00))
  varphi_th_00 <- varphi_th_00 %x% matrix(1,1,nb_strikes)

  array_varphi_th_aa <- array(varphi_th_aa,c(TT,nb_strikes,H))
  array_varphi_th_00 <- array(varphi_th_00,c(TT,nb_strikes,H))

  array_S <- array(S,c(TT,nb_strikes,H))

  matrix_K <- matrix(S,ncol=1) %*% matrix(K_over_S,nrow=1)
  array_K  <- array(matrix_K,c(TT,nb_strikes,H))

  Calls <- array_S * array_varphi_th_aa - array_K * array_varphi_th_00 -
    array_S * G_0a0a + array_K * G_000a

  Puts <- array_K * G_000a - array_S * G_0a0a

  return(
    list(Calls = Calls,
         Puts = Puts)
  )
}



#' Price inflation caps and floors
#'
#' Prices inflation caplets/floorlets and the cumulative caps/floors built from
#' them.
#'
#' @param W Matrix of state-vector values.
#' @param H Maximum maturity in model periods.
#' @param ell Indexation lag.
#' @param all_K Vector of strikes expressed at the model frequency.
#' @param Pi_t_minus_ell Optional lagged price index level.
#' @param psi One-step Laplace transform of the state vector.
#' @param parameterization List containing the model, short-rate parameters, and
#'   inflation loadings.
#' @param max_x,dx_statio,min_dx,nb_x1 Numerical integration settings.
#'
#' @return A list containing `Caplet`, `Floorlet`, `Caps`, `Floors`, and the
#'   maturity grid used in the accumulation.
#'
#' @examples
#' model <- list(n_w = 1, mu = matrix(0, 1, 1), Phi = matrix(0.9, 1, 1),
#'               Sigma12 = matrix(0.1, 1, 1))
#' W <- matrix(c(0, 0.1), 2, 1)
#' res <- price_Inflation_caps_floors(
#'   W = W, H = 2, all_K = c(0.01, 0.02),
#'   psi = psi.GaussianVAR,
#'   parameterization = list(
#'     model = model, xi0 = 0, xi1 = 0,
#'     mu_pi0 = 0.005, mu_pi1 = 0.1
#'   )
#' )
#' dim(res$Caps)
#'
#' @export
price_Inflation_caps_floors <- function(W, # Values of state vector (T x n)
                                        H, # maximum maturity, in model periods
                                        ell = 0, # Indexation lag
                                        all_K, # vector of strikes, expressed at model frequency
                                        # e.g.: quarter. data and K=.01 means annualized strike of 4%
                                        Pi_t_minus_ell = NaN, # Indexation lag multiplicative factor
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
  #     - "mu_pi0" and "mu_pi1" (parameterization of inflation rate, not annualized)

  nb_strikes <- length(all_K)
  n_w <- dim(W)[2]

  xi0 <- parameterization$xi0
  xi1 <- parameterization$xi1

  mu_pi0 <- parameterization$mu_pi0
  mu_pi1 <- parameterization$mu_pi1

  TT <- dim(W)[1] # number of dates.
  vec1TT <- matrix(1,TT,1)

  model <- parameterization$model

  if((ell != 0)&is.na(Pi_t_minus_ell[1])){
    stop("Pi_t_minus_ell cannot be NaN if ell > 0.")
  }

  if((ell == 0) & is.na(Pi_t_minus_ell[1])){
    Pi_t_minus_ell <- matrix(1,TT,1)
  }

  varphi <- function(x,parameterization,H){
    # This function is an argument of truncated.payoff

    u <- parameterization$u
    v <- parameterization$v # determine thresholds (dimension nw x 1)
    i.v.x <- matrix(v,ncol=1) %*% matrix(c(1i*x),nrow=1)

    ell <- parameterization$ell

    nw <- dim(i.v.x)[1]
    q  <- dim(i.v.x)[2]
    U <- matrix(u, nw, q) + i.v.x

    res <- compute_AB_thk(xi0 = parameterization$xi0,
                          xi1 = parameterization$xi1,
                          u = U,
                          H = H, # maximum maturity, in model periods
                          k = ell, # indexation lag of the payoff
                          psi = psi,
                          psi.parameterization = parameterization$model,
                          u2 = NaN,
                          x = 1)
    return(res)
  }

  b.matrix <- matrix(1:H,ncol=1) %*% matrix(all_K, nrow=1) -
    matrix((1:H)-ell,ncol=1) %*% matrix(mu_pi0, 1, nb_strikes) -
    log(Pi_t_minus_ell)

  parameterization$ell <- ell
  parameterization$u   <- NaN

  # Use of "truncated.payoff" for caps:
  parameterization$v <- matrix(mu_pi1,ncol=1)
  parameterization$u <- matrix(mu_pi1,ncol=1)
  res_truncated1 <- truncated.payoff(W,
                                     b.matrix = b.matrix,
                                     varphi = varphi,
                                     parameterization = parameterization,
                                     max_x = max_x,
                                     dx_statio = dx_statio,
                                     min_dx = min_dx,
                                     nb_x1 = nb_x1)
  parameterization$u <- matrix(0*mu_pi1, ncol = 1)
  res_truncated2 <- truncated.payoff(W,
                                     b.matrix = b.matrix,
                                     varphi = varphi,
                                     parameterization = parameterization,
                                     max_x = max_x,
                                     dx_statio = dx_statio,
                                     min_dx = min_dx,
                                     nb_x1 = nb_x1)

  # Computation of inflation swaps:

  res_mu_pi1 <- compute_AB_thk(xi0 = parameterization$xi0,
                               xi1 = parameterization$xi1,
                               u = cbind(mu_pi1,0*mu_pi1),
                               H = H, # maximum maturity, in model periods
                               k = ell, # indexation lag of the payoff
                               psi = psi,
                               psi.parameterization = parameterization$model,
                               u2 = NaN,
                               x = 1)
  # u1 = mu_pi1:
  A1 <- matrix(res_mu_pi1$A[,1,],n_w,H)
  B1 <- matrix(res_mu_pi1$B[,1,],H,1)
  # u1 = 0*mu_pi1:
  A0 <- matrix(res_mu_pi1$A[,2,],n_w,H)
  B0 <- matrix(res_mu_pi1$B[,2,],H,1)

  varphi_th_1 <- exp(W %*% A1 +
                       vec1TT %*% t(B1 + mu_pi0 * matrix((1:H)-ell,H,1)))
  # Along the way, save swap rates:
  swaps <-  log(Pi_t_minus_ell * varphi_th_1) / matrix(1:H,TT,H,byrow = TRUE)
  varphi_th_1 <- varphi_th_1 %x% matrix(1,1,nb_strikes)
  varphi_th_0 <- exp(W %*% A0 + vec1TT %*% t(B0))
  varphi_th_0 <- varphi_th_0 %x% matrix(1,1,nb_strikes)

  array_varphi_th_1 <- array(varphi_th_1,c(TT,nb_strikes,H))
  array_varphi_th_0 <- array(varphi_th_0,c(TT,nb_strikes,H))

  matrix_K <- t(matrix(all_K,nb_strikes,TT))
  array_K  <- array(matrix(1:H,nrow=1) %x% matrix_K,
                    c(TT,nb_strikes,H))
  array_K <- exp(array_K)

  matrix_h_mu_pi0 <- matrix(((1:H)-ell)*mu_pi0,TT,H,byrow=TRUE)
  array_h_mu_pi0 <- array(matrix_h_mu_pi0 %x% matrix(1,1,nb_strikes),
                          c(TT,nb_strikes,H))
  array_h_mu_pi0 <- exp(array_h_mu_pi0)

  array_Pi_t_minus_ell  <- array(Pi_t_minus_ell,c(TT,nb_strikes,H))

  Caps <-
    array_Pi_t_minus_ell * array_varphi_th_1 - array_K * array_varphi_th_0 +
    array_K*res_truncated2 - array_Pi_t_minus_ell*array_h_mu_pi0*res_truncated1

  Floors <- array_K*res_truncated2 -
    array_Pi_t_minus_ell*array_h_mu_pi0*res_truncated1

  dimnames(Caps) <- list(
    time     = paste0("t", 1:TT),
    strike   = paste0(100*all_K,"%"),
    maturity = paste0("matur_",1:H,"_periods")
  )
  dimnames(Floors) <- dimnames(Caps)

  return(
    list(Caps = Caps,
         Floors = Floors,
         swaps = swaps)
  )
}
