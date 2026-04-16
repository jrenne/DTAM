#' Compute first- and second-moment dynamics from an affine transform
#'
#' Numerically differentiates a one-step affine transform to recover the
#' conditional mean and variance of the state vector, together with the implied
#' unconditional moments.
#'
#' @param psi One-step conditional Laplace transform.
#' @param model Model parameter list. It must include `n_w`, the dimension of
#'   the state vector.
#' @param du Finite-difference step used in the numerical derivatives.
#'
#' @return A list containing `mu`, `Phi`, `Theta0`, `Theta1`, `Ew`, and `Vw`.
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
#' res$mu
#' res$Phi
#' res$Ew
#'
#' @export
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
