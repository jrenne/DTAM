solve_PH_TSM <- function(model,
                         max_iter = 100){
  # Solve preferred-habitat-based (PH) term structure model (TSM)
  # It can handle 2 yield curves (_star for the parameters defining the 2nd one)

  # Detect whether one or two yield curves:
  indic_two_curves <- !is.null(model$alpha_star)

  # Maximum maturity:
  H <- length(model$xi)

  # Specification of state vector w_t:
  mu    <- matrix(model$mu,ncol=1)
  Phi   <- model$Phi
  Sigma <- model$Sigma # conditional covariance matrix of w_t
  n     <- sqrt(length(model$Phi)) # number of state variables

  # Net supply specification:
  beta   <- matrix(model$beta,ncol=1)
  alpha  <- matrix(model$alpha,nrow=H)
  xi     <- model$xi
  xi_bar <- matrix(c(xi)/(1:H),ncol=1)

  if(indic_two_curves){# a second yield curve
    beta_star   <- matrix(model$beta_star,ncol=1)
    alpha_star  <- matrix(model$alpha_star,nrow=H)
    xi_star     <- model$xi_star
    xi_bar_star <- matrix(c(xi_star)/(1:H),ncol=1)

    mu_e0 <- model$mu_e0
    mu_e1 <- matrix(model$mu_e1,ncol=1)
  }

  # Specification of short-term rate:
  a1 <- matrix(model$a1,ncol=1)
  b1 <- model$b1

  # Arbitrageurs:
  gamma <- model$gamma

  Gamma <- matrix(0,H,H)
  Gamma[2:H,1:(H-1)] <- diag(H-1)

  vec_1H <- matrix(1,H,1)
  vec_1n <- matrix(1,n,1)

  Id_Gamma_1    <- solve(diag(H) - Gamma)
  Id_PhiGamma_1 <- solve(diag(H*n) - t(Phi) %x% Gamma)

  # Initialization (Sigma = 0):
  A <- matrix(Id_PhiGamma_1 %*% c(-vec_1H %*% t(a1)),H,n)
  Theta <- Gamma %*% A
  B <- Id_Gamma_1 %*% (-b1*vec_1H + Theta %*% mu)
  if(indic_two_curves){# for the second yield curve, if any
    A_star <- matrix(Id_PhiGamma_1 %*% c(-vec_1H %*% t(a1)),H,n)
    Theta_star <- Gamma %*% A_star + vec_1H %*% t(mu_e1)
    B_star <- Id_Gamma_1 %*% (-b1*vec_1H + Theta_star %*% mu)
  }

  for(i in 1:max_iter){

    mathcalB <- beta  + xi_bar * B
    mathcalA <- alpha + (xi_bar %*% t(vec_1n)) * A
    if(indic_two_curves){
      mathcalB_star <- beta_star  + xi_bar_star * B_star
      mathcalA_star <- alpha_star + (xi_bar_star %*% t(vec_1n)) * A_star
    }

    Theta <- Gamma %*% A
    if(indic_two_curves){
      # two curves
      Theta_star <- Gamma %*% A_star + vec_1H %*% t(mu_e1)
      lambda <- gamma * (t(Theta) %*% mathcalB + t(Theta_star) %*% mathcalB_star)
      Lambda <- gamma * (t(Theta) %*% mathcalA + t(Theta_star) %*% mathcalA_star)
    }else{
      # single curve
      lambda <- gamma * t(Theta) %*% mathcalB
      Lambda <- gamma * t(Theta) %*% mathcalA
    }

    muQ  <- mu -  Sigma %*% lambda
    PhiQ <- Phi - Sigma %*% Lambda

    M <- - vec_1H %*% t(a1) -  Theta %*% Sigma %*% Lambda
    A <- matrix(Id_PhiGamma_1 %*% c(M),H,n)
    B <- Id_Gamma_1 %*%
      (-b1*vec_1H + Theta %*% muQ +
         .5 * ((Theta %x% t(vec_1n))*(t(vec_1n) %x% Theta)) %*% c(Sigma) )

    if(indic_two_curves){
      M_star <- - vec_1H %*% t(a1) - Theta_star %*% Sigma %*% Lambda
      A_star <- matrix(Id_PhiGamma_1 %*% c(M_star),H,n)
      B_star <- Id_Gamma_1 %*%
        (-b1*vec_1H + mu_e0*vec_1H + Theta_star %*% muQ +
           .5 * ((Theta_star %x% t(vec_1n))*(t(vec_1n) %x% Theta_star)) %*% c(Sigma) )
    }
  }

  a <- - A / matrix(1:H,H,n)
  b <- - B / matrix(1:H,H,1)

  if(indic_two_curves){
    a_star <- - A_star / matrix(1:H,H,n)
    b_star <- - B_star / matrix(1:H,H,1)
  }else{
    A_star <- NaN
    B_star <- NaN
    a_star <- NaN
    b_star <- NaN
  }

  # Determine specification of net supply for short maturity:
  # (to ensure that sum(z)=1)
  if(indic_two_curves){
    aux <- matrix(mathcalA[2:H,],H-1,n)
    alpha[1,] <- - c(apply(aux,2,sum)) -
      c(apply(mathcalA_star[1:H,],2,sum)) - c(xi[1] * a1)
    beta[1]   <- 1 - sum(mathcalB[2:H]) - sum(mathcalB_star[1:H]) - xi[1] * b1
  } else{
    aux <- matrix(mathcalA[2:H,],H-1,n)
    alpha[1,] <- - c(apply(aux,2,sum)) - c(xi[1] * a1)
    beta[1]   <- 1 - sum(mathcalB[2:H]) - xi[1] * b1
  }

  # Specification of the prices of risk
  lambda0 <- - gamma * t(Theta) %*% mathcalB
  lambda1 <- - gamma * t(Theta) %*% mathcalA

  return(list(A=A, B=B, a=a, b=b,
              alpha=alpha,beta=beta,
              A_star=A_star, B_star=B_star, a_star=a_star, b_star=b_star,
              muQ=muQ, PhiQ=PhiQ,
              lambda0=lambda0, lambda1=lambda1))
}
