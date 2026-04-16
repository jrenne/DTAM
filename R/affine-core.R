#' Reverse multi-horizon Laplace transform recursion
#'
#' Computes the reverse-order multi-horizon Laplace transform recursion for a
#' Markov state process whose one-step conditional Laplace transform is known.
#'
#' Specifically, the function evaluates objects of the form
#' \deqn{
#' \mathbb{E}_t \left[
#'   \exp\left(
#'     u_2^\prime w_{t+1} + \cdots + u_2^\prime w_{t+h-1} + u_1^\prime w_{t+h}
#'   \right)
#' \right].
#' }{
#' E_t[exp(u2' w_{t+1} + ... + u2' w_{t+h-1} + u1' w_{t+h})].
#' }
#'
#' @param psi One-period conditional Laplace transform of the state process.
#' @param u1 Loading matrix applied at the terminal horizon.
#' @param u2 Loading matrix applied at intermediate horizons. If omitted, the
#'   function uses `u2 = u1`.
#' @param H Maximum horizon, in model periods.
#' @param psi.parameterization List of parameters passed to `psi`.
#'
#' @return A list containing arrays `A`, `B`, `a`, and `b`. The arrays `A` and
#'   `B` provide the affine coefficients of the multi-horizon Laplace transform,
#'   while `a` and `b` store the same quantities divided by the horizon.
#'
#' @details
#' The function can process several loading vectors in parallel: if `u1` and
#' `u2` are `n x k`, the recursion is computed simultaneously for the `k`
#' columns.
#'
#' This recursion is a core building block for affine pricing formulas in the
#' package, including bond pricing and several recursive-utility applications.
#'
#' @references
#' Monfort, A., Pegoraro, F., Renne, J.-P., and Roussellet, G. (2026).
#' *Asset Pricing with Discrete-Time Affine Processes*.
#'
#' @examples
#' # Example adapted from the affine-process chapter of the companion Bookdown.
#' model <- list(
#'   mu = matrix(c(0.00, 0.02), ncol = 1),
#'   Phi = matrix(c(0.90, 0.00,
#'                  0.05, 0.85), nrow = 2, byrow = TRUE),
#'   Sigma = diag(c(0.01, 0.02))
#' )
#'
#' U <- matrix(c(0, 1,
#'               -1, 1), nrow = 2)
#'
#' AB <- reverse.MHLT(
#'   psi = psi.GaussianVAR,
#'   u1 = U,
#'   H = 5,
#'   psi.parameterization = model
#' )
#'
#' dim(AB$A)
#' dim(AB$B)
#'
#' @export
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

#' One-step Laplace transform for a Gaussian VAR
#'
#' Computes the one-period conditional Laplace transform associated with a
#' Gaussian VAR(1) state process.
#'
#' The state dynamics are
#' \deqn{
#' w_t = \mu + \Phi w_{t-1} + \varepsilon_t,
#' \qquad \varepsilon_t \sim \mathcal{N}(0, \Sigma).
#' }{
#' w_t = mu + Phi w_{t-1} + eps_t, eps_t ~ N(0, Sigma).
#' }
#'
#' @param u Matrix of Laplace-transform loadings. Each column corresponds to one
#'   transform evaluation.
#' @param psi.parameterization List containing `mu`, `Phi`, and `Sigma`.
#'
#' @return A list with two entries:
#'   `a`, the coefficient multiplying the lagged state in the affine transform,
#'   and `b`, the constant term.
#'
#' @details
#' If `u` has `k` columns, the function computes the `k` transforms in parallel.
#'
#' This function is often used as the `psi` argument of `reverse.MHLT()`,
#' `compute_AB_classical()`, and `compute_AB_thk()` in Gaussian affine models.
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
#' u <- matrix(c(1, 0,
#'               0, 1), nrow = 2)
#' psi.GaussianVAR(u, model)
#'
#' @export
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




