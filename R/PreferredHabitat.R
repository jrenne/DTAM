solve_PH_TSM <- function(model,
                         max_iter = 100,
                         zeta_bar = 3){
  # Solve preferred-habitat-based (PH) term structure model (TSM)
  # It can handle 2 yield curves (_star for the parameters defining the 2nd one)

  # Detect whether one or two yield curves:
  indic_two_curves <- !is.null(model$alpha_star)

  # Detect whether stock prices:
  indic_stocks <- !is.null(model$omega)

  if(indic_stocks&indic_two_curves){
    print("ERROR: The code cannot handle 2 yield curves & stocks simultaneously.")
    return(0)
  }

  # Maximum maturity:
  H <- length(model$xi)

  # Specification of state vector w_t:
  mu    <- matrix(model$mu,ncol=1)
  Phi   <- model$Phi
  Sigma <- model$Sigma # conditional covariance matrix of w_t
  n     <- sqrt(length(model$Phi)) # number of state variables
  Ew    <- solve(diag(n) - Phi) %*% mu # uncondi. expectation

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

    if(!is.null(model$mu_e0)){
      # In that case, mu_e0 and mu_e1 are known, a1_star and b1_star are deduced
      mu_e0 <- model$mu_e0
      mu_e1 <- matrix(model$mu_e1,ncol=1)
    }else{
      # In that case, a1_star and b1_star are known, mu_e0 and mu_e1 are deduced
      a1_star <- model$a1_star
      b1_star <- model$b1_star
    }
  }

  if(indic_stocks){# stocks (or perpetuity) in the model
    omega <- model$omega
    mu_K0 <- model$mu_K0
    mu_K1 <- model$mu_K1

    beta_K   <- model$beta_K
    alpha_K  <- matrix(model$alpha_K,ncol=1)
    xi_K     <- model$xi_K
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
    if(is.null(model$mu_e0)){
      mu_e1 <- solve(t(Phi)) %*% (a1 - a1_star)
      mu_e0 <- c(b1 - b1_star - t(mu_e1) %*% mu)
    }
    # A_star <- matrix(Id_PhiGamma_1 %*% c(
    #   -vec_1H %*% t(a1) + (vec_1H %*% t(mu_e1)) %*% Phi
    # ),H,n)
    # Theta_star <- Gamma %*% A_star + vec_1H %*% t(mu_e1)
    # B_star <- Id_Gamma_1 %*% (-b1*vec_1H + Theta_star %*% mu)

    A_star <- 0*A
    B_star <- 0*B
    Theta_star <- 0*Theta
  }

  if(indic_stocks){
    kappa1 <- exp(zeta_bar)/(1+exp(zeta_bar))
    kappa0 <- log(1+exp(zeta_bar)) - kappa1 * zeta_bar
    mu_zeta1 <- solve((1 + omega*xi_K)*diag(n) - kappa1*t(Phi)) %*%
      (- a1 - omega*alpha_K + t(Phi)%*%mu_K1)
    mu_zeta0 <- 1/(1 - kappa1 + omega*xi_K)*
      (- b1 - omega*beta_K + kappa0 + mu_K0 + t(kappa1*mu_zeta1 + mu_K1)%*%mu)
    zeta_bar <- c(mu_zeta0 + t(mu_zeta1) %*% Ew)
  }

  for(i in 1:max_iter){

    mathcalB <- beta  + xi_bar * B
    mathcalA <- alpha + (xi_bar %*% t(vec_1n)) * A
    if(indic_two_curves){
      mathcalB_star <- beta_star  + xi_bar_star * B_star
      mathcalA_star <- alpha_star + (xi_bar_star %*% t(vec_1n)) * A_star
    }
    if(indic_stocks){
      mathcalB_K <- beta_K  + xi_K * mu_zeta0
      mathcalA_K <- matrix(alpha_K + xi_K * mu_zeta1,nrow=1)
    }

    # Compute Thetas and prices of risk:
    Theta <- Gamma %*% A
    if(indic_two_curves){
      # two curves
      Theta_star <- Gamma %*% A_star + vec_1H %*% t(mu_e1)
      lambda0 <- - gamma * (t(Theta) %*% mathcalB + t(Theta_star) %*% mathcalB_star)
      lambda1 <- - gamma * (t(Theta) %*% mathcalA + t(Theta_star) %*% mathcalA_star)
      if(is.null(model$mu_e0)){
        mu_e1 <- solve(t(Phi) + t(lambda1) %*% Sigma) %*% (a1 - a1_star)
        mu_e0 <- c(b1 - b1_star - t(mu_e1) %*% (mu + Sigma %*% lambda0) -
                     .5 * t(mu_e1) %*% Sigma %*% mu_e1)
        Theta_star <- Gamma %*% A_star + vec_1H %*% t(mu_e1)
      }
    }else if(indic_stocks){
      # stock returns
      Theta_K <- matrix(kappa1*mu_zeta1 + mu_K1,nrow=1)
      lambda0 <- - gamma * (t(Theta) %*% mathcalB + t(Theta_K) %*% mathcalB_K)
      lambda1 <- - gamma * (t(Theta) %*% mathcalA + t(Theta_K) %*% mathcalA_K)
    }else{
      # single curve
      lambda0 <- - gamma * t(Theta) %*% mathcalB
      lambda1 <- - gamma * t(Theta) %*% mathcalA
    }

    muQ  <- mu + Sigma %*% lambda0
    PhiQ <- Phi + Sigma %*% lambda1

    M <- - vec_1H %*% t(a1) +  Theta %*% Sigma %*% lambda1
    A <- matrix(Id_PhiGamma_1 %*% c(M),H,n)
    B <- Id_Gamma_1 %*%
      (-b1*vec_1H + Theta %*% muQ +
         .5 * ((Theta %x% t(vec_1n))*(t(vec_1n) %x% Theta)) %*% c(Sigma) )

    if(indic_two_curves){
      M_star <- - vec_1H %*% t(a1) + Theta_star %*% Sigma %*% lambda1 +
        (vec_1H %*% t(mu_e1)) %*% Phi
      A_star <- matrix(Id_PhiGamma_1 %*% c(M_star),H,n)
      B_star <- Id_Gamma_1 %*%
        (-b1*vec_1H + mu_e0*vec_1H + Theta_star %*% muQ +
           .5 * ((Theta_star %x% t(vec_1n))*(t(vec_1n) %x% Theta_star)) %*% c(Sigma) )
    }
    if(indic_stocks){
      kappa1 <- exp(zeta_bar)/(1+exp(zeta_bar))
      kappa0 <- log(1+exp(zeta_bar)) - kappa1 * zeta_bar
      mu_zeta1 <- solve((1 + omega*xi_K)*diag(n) - kappa1*t(PhiQ)) %*%
        (- a1 - omega*alpha_K + t(PhiQ)%*%mu_K1)
      mu_zeta0 <- 1/(1 - kappa1 + omega*xi_K)*
        (- b1 - omega*beta_K + kappa0 + mu_K0 +
           1/2 * t(kappa1*mu_zeta1 + mu_K1) %*% Sigma %*% t(kappa1*mu_zeta1 + mu_K1) +
           t(kappa1*mu_zeta1 + mu_K1)%*%muQ)
      zeta_bar <- c(mu_zeta0 + t(mu_zeta1) %*% Ew)
      #print(zeta_bar)
    }
  }

  a <- - A / matrix(1:H,H,n)
  b <- - B / matrix(1:H,H,1)

  if(indic_two_curves){
    a_star <- - A_star / matrix(1:H,H,n)
    b_star <- - B_star / matrix(1:H,H,1)
    if(!is.null(model$mu_e0)){
      a1_star <- matrix(a_star[1,],ncol=1)
      b1_star <- b_star[1]
    }
  }else{
    A_star  <- NaN
    B_star  <- NaN
    a_star  <- NaN
    b_star  <- NaN
    a1_star <- NaN
    b1_star <- NaN
    mu_e0   <- NaN
    mu_e1   <- NaN
  }

  if(indic_stocks){
    C_K <- kappa0 + (kappa1-1)*mu_zeta0 + mu_K0
    D_K <- - t(mu_zeta1)
    Theta_K <- t(kappa1*mu_zeta1 + mu_K1)
  }else{
    mu_zeta0 <- NaN
    mu_zeta1 <- NaN
    kappa0 <- NaN
    kappa1 <- NaN

    C_K     <- NaN
    D_K     <- NaN
    Theta_K <- NaN
  }

  # Determine specification of net supply for short maturity:
  # (to ensure that sum(z)=1)
  aux <- matrix(mathcalA[2:H,],H-1,n)
  if(indic_two_curves){
    alpha[1,] <- - c(apply(aux,2,sum)) -
      c(apply(mathcalA_star[1:H,],2,sum)) - c(xi[1] * a1)
    beta[1]   <- 1 - sum(mathcalB[2:H]) - sum(mathcalB_star[1:H]) - xi[1] * b1
  }else if(indic_stocks){
    alpha[1,] <- - c(apply(aux,2,sum)) - c(xi[1] * a1) - c(mathcalA_K)
    beta[1]   <- 1 - sum(mathcalB[2:H]) - xi[1] * b1 - mathcalB_K
  }else{
    alpha[1,] <- - c(apply(aux,2,sum)) - c(xi[1] * a1)
    beta[1]   <- 1 - sum(mathcalB[2:H]) - xi[1] * b1
  }

  return(list(A=A, B=B, a=a, b=b,
              alpha=alpha,beta=beta,
              A_star=A_star, B_star=B_star, a_star=a_star, b_star=b_star,
              muQ=muQ, PhiQ=PhiQ,
              lambda0=lambda0, lambda1=lambda1, zeta_bar=zeta_bar,
              mu_zeta0=mu_zeta0, mu_zeta1=mu_zeta1, C_K=C_K, D_K=D_K,
              Theta_K=Theta_K,
              a1_star=a1_star,b1_star=b1_star,mu_e0=mu_e0,mu_e1=mu_e1))
}
