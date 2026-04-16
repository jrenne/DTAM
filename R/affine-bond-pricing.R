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
#' @return A list with arrays `A`, `B`, `a`, and `b`.
#'
#' @examples
#' model <- list(
#'   mu = matrix(c(0.01, 0.02), ncol = 1),
#'   Phi = matrix(c(0.9, 0.05,
#'                  0.00, 0.85), nrow = 2, byrow = TRUE),
#'   Sigma = diag(c(0.02, 0.01))
#' )
#' res <- compute_AB_classical(
#'   xi0 = 0.01,
#'   xi1 = matrix(c(0.1, 0.2), ncol = 1),
#'   H = 4,
#'   psi = psi.GaussianVAR,
#'   psi.parameterization = model
#' )
#' dim(res$A)
#'
#' @export
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

#' Compute affine coefficients for delayed-horizon payoffs
#'
#' Computes affine pricing coefficients for payoffs involving the state vector
#' at a delayed horizon and, optionally, sums of intermediate future states.
#'
#' The function is used for claims of the form
#' \deqn{
#' \exp\left(
#'   u^\prime \left(w_{t+h-k} + x[w_{t+h-k-1} + \cdots + w_{t+1}]\right)
#'   + u_2^\prime (w_{t+h} + \cdots + w_{t+1})
#' \right),
#' }{
#' exp(u' (w_{t+h-k} + x[w_{t+h-k-1} + ... + w_{t+1}]) +
#'   u2' (w_{t+h} + ... + w_{t+1})).
#' }
#'
#' @param xi0 Scalar intercept in the short-rate equation.
#' @param xi1 Loading vector in the short-rate equation.
#' @param u Loading matrix applied to the delayed-horizon state.
#' @param H Maximum maturity, in model periods.
#' @param k Indexation lag of the payoff, in model periods.
#' @param psi Conditional Laplace transform of the state vector.
#' @param psi.parameterization List of parameters passed to `psi`.
#' @param u2 Optional loading matrix applied to the cumulative future state
#'   terms. If left as `NaN`, the function sets `u2 = 0`.
#' @param x Scalar multiplier applied to the intermediate-state sum.
#'
#' @return A list with arrays `A` and `B`, the affine coefficients of the
#'   date-`t` price.
#'
#' @details
#' The function returns `NaN` coefficients for horizons `h < k`, since the
#' delayed payoff is not yet defined at those maturities.
#'
#' This helper is used in option-pricing and indexed-payoff applications in the
#' Bookdown companion project.
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
#'   Sigma = diag(c(0.02, 0.01))
#' )
#'
#' xi0 <- 0
#' xi1 <- matrix(c(1, 0), ncol = 1)
#' u <- matrix(c(0, 1), ncol = 1)
#'
#' res <- compute_AB_thk(
#'   xi0 = xi0,
#'   xi1 = xi1,
#'   u = u,
#'   H = 6,
#'   k = 2,
#'   psi = psi.GaussianVAR,
#'   psi.parameterization = model
#' )
#'
#' dim(res$A)
#' dim(res$B)
#'
#' @export
compute_AB_thk <- function(xi0, xi1,
                           u,
                           H, # maximum maturity,
                           k, # indexation lag of the payoff
                           psi, psi.parameterization,
                           u2 = NaN, x = 0){
  # This function is used to price options.
  # It employs recursive formulas to price payoffs of the type:
  #    exp{u'(w_{t+h-k} + x*[w_{t+h-k-1} + ... + w_{t+1}])
  #        + u2'(w_{t+h}+...+w_{t+1})},
  #    settled on date t+h.
  #    (if u2 = NaN, then u2 <- 0)
  # The risk-free short-term rate is r_t = xi0 + xi1'w_t,
  #    where w_t is the state vector, of dimension n x 1.
  # The function returns a list containing A, B, a, and b.
  # Bond prices are such that the date-t payoff is = exp(B_h + A_h'w_t), where
  #    B_h is a scalar and A_h is a (n x 1) vector.
  # u can be of dimension n x q with q > 1 (n is the dimension of the state)
  # psi and psi.parameterization defined the risk-neutral dynamics of w_t.
  # The function returns NaN for horizons h < k.

  n <- length(xi1) # dimension of state vector
  u1 <- matrix(u,nrow=n) # make sure u has the right dimension

  if(is.na(u2[1])){
    u2 <- 0*u1
  }

  A.h_1 <- (k==0) * u1
  B.h_1 <- 0

  n <- dim(u)[1]
  q <- dim(u)[2]

  A <- array(NaN, c(n, q, H))
  B <- array(NaN, c(1, q, H))
  for (h in 1:H) {
    X <- matrix(-xi1, n, q) * (h > 1) + A.h_1 + u2 + u1 * x * (h > k)
    psi.u <- psi(X, psi.parameterization)
    A.h <- matrix(psi.u$a, n, q) + u1 * (h == k)
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

#' Price interest-rate caps and floors
#'
#' Prices caplets, floorlets, caps, and floors in an affine term-structure model
#' using Fourier inversion and affine pricing recursions.
#'
#' @param W Matrix of state values, with one row per date.
#' @param H Maximum maturity, in model periods.
#' @param tau Maturity of the reference rate, in model periods.
#' @param freq Number of model periods per year.
#' @param all_K Vector of annualized strike rates.
#' @param psi Conditional Laplace transform of the state vector.
#' @param parameterization List containing at least `model`, `xi0`, and `xi1`.
#'   The embedded `model` object must include `n_w`.
#' @param max_x Maximum truncation point of the numerical integration grid.
#' @param dx_statio Constant step size used on the upper tail of the integration
#'   grid.
#' @param min_dx Smallest step size used at the beginning of the grid.
#' @param nb_x1 Number of log-spaced grid points used before switching to
#'   `dx_statio`.
#'
#' @return A list containing:
#'   `A_tau`, `B_tau`, `a_tau`, `b_tau` for the affine representation of the
#'   reference rate, plus arrays `Caplet`, `Floorlet`, `Caps`, and `Floors`.
#'
#' @details
#' The function first computes the affine representation of the reference rate,
#' then prices caplets and floorlets through repeated calls to
#' `truncated.payoff()`, and finally aggregates them into caps and floors across
#' maturities.
#'
#' The entries of `all_K` are annualized strikes, whereas the underlying model
#' is expressed at the model frequency set by `freq`.
#'
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#'
#' @examples
#' model <- list(
#'   alpha = matrix(c(0.1, 0.2), ncol = 1),
#'   nu = matrix(c(0.2, 0.3), ncol = 1),
#'   mu = matrix(c(1, 1), ncol = 1),
#'   beta = matrix(c(0.9, 0,
#'                   0.1, 0.8), 2, 2, byrow = TRUE),
#'   n_w = 2
#' )
#'
#' set.seed(123)
#' W <- simul_VARG(model, nb_periods = 30)
#' W <- t(W)
#'
#' parameterization <- list(
#'   model = model,
#'   xi0 = 0.01 / 12,
#'   xi1 = matrix(c(0.5, 0.2), ncol = 1)
#' )
#'
#' res <- price_IR_caps_floors(
#'   W = W,
#'   H = 12,
#'   tau = 3,
#'   freq = 12,
#'   all_K = c(0.01, 0.02),
#'   psi = psi.VARG,
#'   parameterization = parameterization,
#'   max_x = 500,
#'   dx_statio = 2,
#'   min_dx = 1e-4,
#'   nb_x1 = 300
#' )
#'
#' names(res)
#' dim(res$Caps)
#'
#' @export
