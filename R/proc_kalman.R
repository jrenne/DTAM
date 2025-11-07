
Kalman_filter <- function(Y_t,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,
                          indic_pos=0,
                          Rfunction=Rf, Qfunction=Qf,
                          reconciliationf = function(x){x} # In case
){
  # Measurement equations:
  #    y_t   = mu_t + G * rho_t + M * eps_t
  # Transition equations:
  #    rho_t = nu_t + H * rho_t-1 + N * xi_t
  # Y_t is a T*ny-dimensional vector
  # nu_t and mu_t depend on past only (predetermined at date t)
  #
  # The function returns:
  # .$r   filtered variables (rho_t|t) size T*nr
  # .$S  SIGMA t|t T*(nr*nr)

  # Notes:
  # 'indic_pos' is a vector of binary variables indicating if
  #     latent factors have to be >=0
  # If 'reconciliationf' is a function, then 'rho_tt' may be modified,
  #     in a way defined by this function. That is used in the QKF, to ensure
  #     that there is consistency in the augmented state vector (x,vec(xx')).

  # Number of observed variables:
  ny = NCOL(Y_t)
  # Number of unobs. variables:
  nr = NCOL(G)
  # Number of time periods:
  T = NROW(Y_t)

  # loglik.vector will contain the vector of date-specific log-likelihood
  loglik.vector <- NULL

  #Initilize output matrices:
  rho_tt    = matrix(0,T,nr)
  rho_tp1_t = matrix(0,T,nr)
  y_tp1_t   = matrix(0,T,ny)

  Sigma_tt    = matrix(0,T,nr*nr)
  Sigma_tp1_t = matrix(0,T,nr*nr)
  Omega_tt    = matrix(0,T,ny*ny)
  Omega_tp1_t = matrix(0,T,ny*ny)

  #Initilize log-Likelihood:
  logl = ny*T/2*log(2*pi)

  #Kalman algorithm:

  for (t in 1:T){

    # print(t)

    # ==========================================================================
    # Forecasting step (between t-1 and t):
    # ==========================================================================
    if(t==1){
      rho_tp1_t[1,] = nu_t[1,] + t(H %*% rho_0)
      R = Rfunction(M,rho_0) #Rfunction defined below (measurement eq.)
      Q = Qfunction(N,rho_0) #Qfunction defined below (transition eq.)
    } else {
      rho_tp1_t[t,] = nu_t[t,] + t(H %*% rho_tt[t-1,])
      R = Rfunction(M,rho_tt[t-1,],t)
      Q = Qfunction(N,rho_tt[t-1,],t)
    }
    if(sum(indic_pos==1)>0){ # Some latent variables imposed to be >= 0
      rho_tp1_t[t,indic_pos==1] = pmax(rho_tp1_t[t,indic_pos==1],0)
    }
    y_tp1_t[t,] = mu_t[t,] + t(G %*% rho_tp1_t[t,])

    if(t==1){
      aux_Sigma_tp1_t = Q + H %*% Sigma_0 %*% t(H)
    } else {
      aux_Sigma_tp1_t = Q + H %*% matrix(Sigma_tt[t-1,],nr,nr) %*% t(H)
    }

    Sigma_tp1_t[t,] = matrix( aux_Sigma_tp1_t ,1,nr*nr)
    omega           = R + G %*% aux_Sigma_tp1_t %*% t(G)
    Omega_tp1_t[t,] = matrix(omega,1,ny*ny)


    # ==========================================================================
    # Updating step
    # ==========================================================================

    # Detect observed variables:
    vec.obs.indices <- which(!is.na(Y_t[t,])) # indices of observed variables
    ny.aux <- length(c(vec.obs.indices)) # number of observed variables
    # Resize matrices accordingly:
    G.aux <- matrix(G[vec.obs.indices,],nrow=ny.aux)
    R.aux <- R[vec.obs.indices,]
    omega <- omega[vec.obs.indices,]
    if(class(R.aux)[1]=="numeric"){
      R.aux <- R.aux[vec.obs.indices]
      omega <- omega[vec.obs.indices]
    }else{
      R.aux <- R.aux[,vec.obs.indices]
      omega <- omega[,vec.obs.indices]
    }
    R.aux <- matrix(R.aux,ny.aux,ny.aux)

    #Compute gain K:
    if(dim(G.aux)[1]>0){
      K = aux_Sigma_tp1_t %*% t(G.aux) %*% MASS::ginv(R.aux + G.aux %*% aux_Sigma_tp1_t %*% t(G.aux))

      lambda_t   = Y_t[t,] - y_tp1_t[t,]
      lambda_t <- matrix(lambda_t[vec.obs.indices],ncol=1)
      rho_tt[t,] = t( rho_tp1_t[t,] + K %*% lambda_t )
      if(sum(indic_pos==1)>0){
        rho_tt[t,indic_pos==1] = pmax(rho_tt[t,indic_pos==1],0)
      }
      Id         = diag(1,nrow=nr,ncol=nr)
      Sigma_tt[t,] = matrix( (Id - K %*% G.aux) %*% aux_Sigma_tp1_t,1,nr*nr)

      # Computation of log-likelihood:
      if(length(c(omega))==1){
        det.omega <- omega
      }else{
        det.omega <- det(omega)
      }

      # Potential reconciliation between components of rho_tt (QKF):
      y2Bfitted <- matrix(Y_t[t,],ncol=1)
      constant  <- matrix(mu_t[t,],ncol=1)
      rho_tt[t,] <- reconciliationf(rho_tt[t,],
                                    y2Bfitted,constant,G,M)

      loglik.vector <- rbind(loglik.vector,
                             ny.aux/2*log(2*pi) - 1/2*(log(det.omega) +
                                                         t(lambda_t) %*% MASS::ginv(omega) %*% lambda_t)
      )
      logl <- logl + loglik.vector[t]

    }else{
      rho_tt[t,] = t( rho_tp1_t[t,])
      Sigma_tt[t,] = matrix( aux_Sigma_tp1_t,1,nr*nr)
    }

  }

  fitted.obs <- mu_t + rho_tt %*% t(G) # fitted observables

  output = list(r=rho_tt,Sigma_tt=Sigma_tt,loglik=logl,y_tp1_t=y_tp1_t,
                S_tp1_t=Sigma_tp1_t,r_tp1_t=rho_tp1_t,
                loglik.vector = loglik.vector,Omega_tp1_t=Omega_tp1_t,M=M,
                fitted.obs=fitted.obs)
  return(output)
}

Rf <- function(M,RHO,t=0){
  return(M %*% t(M))
}
Qf <- function(N,RHO,t=0){
  return(N %*% t(N))
}

Kalman_smoother <- function(Y_t,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,indic_pos=0,
                            Rfunction=Rf, Qfunction=Qf){

  output_filter <- Kalman_filter(Y_t,nu_t,H,N,mu_t,G,M,Sigma_0,rho_0,indic_pos,
                                 Rfunction,Qfunction)
  rho_tt        <- output_filter$r
  Sigma_tt      <- output_filter$Sigma_tt
  Sigma_tp1_t   <- output_filter$S_tp1_t
  rho_tp1_t     <- output_filter$r_tp1_t
  Omega_tp1_t   <- output_filter$Omega_tp1_t

  # Number of observed variables:
  ny = NCOL(Y_t)
  # Number of unobs. variables:
  nr = NCOL(G)
  # Number of time periods:
  TT = NROW(Y_t)

  if(class(M)[1]=="list"){
    M.aux <- M$sigmas
  }else{
    M.aux <- M
  }

  #R = M%*%t(M)
  #Q = N%*%t(N)

  #Initilize output matrices:
  rho_tT       <- matrix(0,TT,nr)
  rho_tT[TT,]   <- rho_tt[TT,]
  Sigma_tT     <- matrix(0,TT,nr*nr)
  Omega_tT     <- matrix(0,TT,ny*ny) # This is the covariance matrix of Y_t if Y_t is missing
  Sigma_tT[TT,] <- Sigma_tt[TT,]
  Sigma_tT.aux <- matrix(Sigma_tt[TT,],nr,nr)
  R = Rfunction(M,rho_0) #Rfunction defined below (measurement eq.)
  Omega_tT[TT,] <- c(
    R + G %*% Sigma_tT.aux %*% t(G)
  )
  Var.model.implied.obs <- NaN * Y_t # this matrix will contain the variances of G*rho_t
  Var.model.implied.obs[TT,] <- diag(
    G %*% matrix(Sigma_tT[TT,],nr,nr) %*% t(G))

  for(t in seq(TT-1,1,by=-1)){
    F <- matrix(Sigma_tt[t,],nr,nr) %*% t(H) %*% solve(matrix(Sigma_tp1_t[t+1,],nr,nr))
    rho_tT[t,] <- rho_tt[t,] + t( F %*% (rho_tT[t+1,]-rho_tp1_t[t+1,]) )
    if(sum(indic_pos==1)>0){
      rho_tT[t,indic_pos==1] = pmax(rho_tT[t,indic_pos==1],0)
    }
    Sigma_tT.aux <- Sigma_tt[t,] + F %*%
      matrix(Sigma_tT[t+1,]-Sigma_tp1_t[t+1,],nr,nr) %*% t(F)
    Sigma_tT[t,] <- c(Sigma_tT.aux)
    R = Rfunction(M,rho_0) #Rfunction defined below (measurement eq.)
    Omega_tT[t,] <- c(
      R + G %*% Sigma_tT.aux %*% t(G)
    )
  }
  fitted.obs <- mu_t + rho_tT %*% t(G) # fitted observables
  output = list(r_smooth=rho_tT,S_smooth=Sigma_tT,
                loglik=output_filter$loglik,
                loglik.vector=output_filter$loglik.vector,
                Omega_tp1_t=Omega_tp1_t,
                M=M,
                Sigma_tt=Sigma_tt,
                fitted.obs=fitted.obs,G=G,mu_t=mu_t,
                Omega_tT = Omega_tT)
  return(output)
}

QKF <- function(Y_t,QStateSpace,indic_reconciliation=TRUE){
  # The state-space is: --------------------------------------------------------
  # --- Measurement equations (y_t of dimension m x 1):
  # y_t = A + B·x_t + Sum_i{e_i·[x_t'·C_i·x_t]} + eta_t,
  #   with eta_t ~ N(0,Omega), with Omega = M·M'.
  # --- Transition equation (x_t of dimension n x 1):
  # x_t = mu + Phi·x_{t-1} + epsilon_t,
  #   with epsilon_t ~ N(0,Sigma)
  # Remarks: -------------------------------------------------------------------
  # C is an array of dimension n x n x m; the i^th layer is C_i.
  #    The layers C_i are not necessarily symmetric.
  # The dimension of the augmented state vector is n + n*(n+1)/2,
  #    that is, the augmented state vector is (x_t',vech(x_t·x_t')')'
  # ----------------------------------------------------------------------------

  A     <- QStateSpace$A
  B     <- QStateSpace$B
  C     <- QStateSpace$C
  mu    <- QStateSpace$mu
  Phi   <- QStateSpace$Phi
  Sigma <- QStateSpace$Sigma
  M     <- QStateSpace$M

  B <- as.matrix(B)
  if(class(C) != "array"){
    print("Error: C should be an array")
    return(0)
  }

  # Dimensions:
  T <- dim(Y_t)[1]
  m <- dim(B)[1]
  n <- dim(B)[2]

  model <- list(mu=mu,Phi=Phi,Sigma=Sigma)
  # Compute mu and Psi:
  aux <- make_matrices_cond_mean_variance_Quadratic(model,
                                                    indic_compute_V=TRUE)

  # Convert inputs into those required by 'Kalman_filter':
  nu_t <- t(matrix(c(aux$mu_tilde_red),n + n*(n+1)/2,T))
  H    <- aux$Phi_tilde_red
  N    <- list(n=n,
               nu_red = aux$nu_red,
               Psi_red = aux$Psi_red)
  mu_t <- t(matrix(c(A),m,T))

  Cmatrix     <- t(apply(C, 3, as.vector)) # number of columns: n^2
  Kn          <- make_Knx(n) # Kn used to build Cmatrix_red, with n*(n+1)/2 columns.
  Cmatrix_red <- Cmatrix %*% t(Kn) # number of columns: n*(n+1)/2
  G           <- cbind(B,Cmatrix_red)

  M <- as.matrix(M)

  Ex <- matrix(aux$E_red[1:n],ncol=1)
  rho_0   <- c(Ex,(Ex %*% t(Ex))[!upper.tri(Ex %*% t(Ex))])
  Sigma_0 <- aux$V_red

  if(indic_reconciliation){
    resKF <- Kalman_filter(Y_t,nu_t,H,N,mu_t,G,M,
                           Sigma_0=Sigma_0,rho_0=rho_0,
                           indic_pos=0,
                           Rfunction=Rf, Qfunction=Q_QKF,
                           reconciliationf = reconciliationf_QKF)
  }else{
    resKF <- Kalman_filter(Y_t,nu_t,H,N,mu_t,G,M,
                           Sigma_0=Sigma_0,rho_0=rho_0,
                           indic_pos=0,
                           Rfunction=Rf, Qfunction=Q_QKF)
  }

  return(resKF)
}

Q_QKF <- function(N,RHO,t=0){
  nu_red  <- N$nu_red
  Psi_red <- N$Psi_red
  n <- N$n
  Q <- matrix(nu_red + Psi_red %*% RHO,n+n*(n+1)/2,n+n*(n+1)/2)
  return(Q)
}

reconciliationf_QKF <- function(rho_ini,y2Bfitted,constant,G,M){
  Omega_1 <- solve(M %*% t(M))

  # Determine n:
  r <- dim(G)[2]
  D <- 9 + 8*r
  n <- (-3+sqrt(D))/2

  # Linear projection:
  x <- matrix(rho_ini[1:n],ncol=1)
  S <- x %*% t(x)
  w1 <- matrix(c(x,S[!upper.tri(S)]),ncol=1)
  lambda <- y2Bfitted - (constant + G %*% w1)
  ell1 <- c(t(lambda) %*% Omega_1 %*% lambda)

  # Quadratic projection:
  vechS <- rho_ini[(n+1):(n+n*(n+1)/2)]
  vecS <- t(make_Knx(n)) %*% vechS
  S <- matrix(vecS,n,n)
  eigenS <- eigen(S)
  v1 <- matrix(eigenS$vectors[,1],ncol=1)
  lambda1 <- eigenS$values[1]
  x <- sqrt(lambda1) * v1
  S <- x %*% t(x)
  w2 <- matrix(c(x,S[!upper.tri(S)]),ncol=1)
  lambda <- y2Bfitted - (constant + G %*% w2)
  ell2 <- c(t(lambda) %*% Omega_1 %*% lambda)
  if(ell1 < ell2){# better fit with linear projection
    rho_tt <- w1
  }else{# better fit with quadratic projection
    rho_tt <- w2
  }
  return(rho_tt)
}

