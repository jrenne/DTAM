#' VARG transforms, moments, and simulation helpers
#'
#' These functions cover the vector autoregressive gamma (VARG) model used in
#' the credit chapters. They provide the one-step Laplace transform under the
#' physical or transformed measure, unconditional moments, simulation, default
#' probabilities, and CDS pricing quantities.
#'
#' @param u Matrix of transform arguments.
#' @param psi.parameterization,model Lists of model parameters.
#' @param H Maximum horizon.
#' @param indic_delta1 Optional index locating the first default-intensity
#'   component within the state vector.
#' @param X Matrix of state-vector values where the pricing objects are
#'   evaluated.
#' @param nb_periods Number of simulated dates.
#' @param w0 Optional initial state.
#'
#' @return Depending on the function: affine-transform coefficients (`a`, `b`),
#'   unconditional moments, simulated paths, transformed-measure parameters, or
#'   pricing objects such as default probabilities and CDS spreads.
#'
#' @examples
#' model <- list(
#'   alpha = c(0.5, 0.3),
#'   nu = c(1, 1.2),
#'   mu = c(0.2, 0.15),
#'   beta = matrix(c(0.4, 0.1,
#'                   0.0, 0.3), 2, 2, byrow = TRUE),
#'   alpha_w = c(0.1, 0.05),
#'   mu_R0 = 0,
#'   mu_R1 = matrix(0, 2, 1),
#'   xi0 = 0,
#'   xi1 = matrix(0, 2, 1)
#' )
#' res_psi <- psi.VARG(matrix(c(0.05, 0.02), 2, 1), model)
#' dim(res_psi$a)
#' compute_uncondmean_VARG(model)
#'
#' @name VARG_tools
#' @export
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


#' @rdname VARG_tools
#' @export
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


#' @rdname VARG_tools
#' @export
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

#' @rdname VARG_tools
#' @export
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

#' @rdname VARG_tools
#' @export
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

#' @rdname VARG_tools
#' @export
make_VARG_Q <- function(model){
  modelQ <- model
  modelQ$mu    <- model$mu / (1 - model$alpha_w * model$mu)
  modelQ$alpha <- model$alpha * modelQ$mu/model$mu
  modelQ$beta  <- diag(c(modelQ$mu/model$mu)) %*% model$beta
  return(modelQ)
}

#' @rdname VARG_tools
#' @export
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
#' Top-down credit-model transforms and simulation
#'
#' The `psi.*TopDown` functions evaluate the one-step Laplace transform of the
#' latent loss/intensity state used in the top-down credit model. `simul.TopDown`
#' simulates the state vector under the physical measure.
#'
#' @param u,u.y Matrix of transform arguments.
#' @param model List of top-down model parameters.
#' @param psi.y Optional function used to evaluate the transform of the `y`
#'   block.
#' @param nb_periods Number of simulated dates.
#' @param W0 Optional initial state.
#'
#' @return The `psi.*` functions return affine-transform coefficients. The
#'   simulator returns a list with unconditional moments (`EV`) and simulated
#'   paths for `y`, `n`, and the stacked state `w`.
#'
#' @examples
#' model <- list(
#'   alpha.n = c(0.2),
#'   beta.ny = matrix(0.1, 1, 1),
#'   beta.nn = matrix(0.3, 1, 1),
#'   beta.yy = matrix(0.2, 1, 1),
#'   beta.yn = matrix(0.05, 1, 1),
#'   alpha.y = c(0.1),
#'   nu.y = c(1),
#'   mu.y = c(0.4),
#'   alpha.m = matrix(0, 2, 1),
#'   n_w = 2
#' )
#' psi.y.TopDown(matrix(0.05, 1, 1), model)$b
#' sim <- simul.TopDown(model, nb_periods = 10)
#' dim(sim$all_w)
#'
#' @name TopDown_tools
#' @export
