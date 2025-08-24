
make_recessions <- function(col="#AA55AA44",indic_US=TRUE,
                            MIN=-10000,MAX=+10000){
  if(indic_US){
    nber_dates <- nberDates()
    start_recession_dates <- as.Date(as.character(nber_dates[,1]),"%Y%m%d")
    end_recession_dates   <- as.Date(as.character(nber_dates[,2]),"%Y%m%d")
  }else{
    # Euro zone
    start_recession_dates <- as.Date(c("2008-02-01","2011-08-01","2019-11-01"))
    end_recession_dates   <- as.Date(c("2009-05-01","2013-02-01","2020-08-01"))
  }

  for(i in 1:length(start_recession_dates)){
    polygon(c(start_recession_dates[i],start_recession_dates[i],
              end_recession_dates[i],end_recession_dates[i]),
            c(MIN,MAX,MAX,MIN),border=NaN,col=col)
  }
  return(NULL)
}

autocov <- function(X,n){
  X <- as.matrix(X)
  if(dim(X)[1]<dim(X)[2]){
    X <- t(X)
  }
  T <- dim(X)[1]
  X <- X - matrix(1,T,1) %*% matrix(apply(X,2,mean),nrow=1)
  X.1 <- X[1:(T-n),]
  X.2 <- X[(n+1):T,]
  return(1/T * t(X.1) %*% X.2)
}

NW.LongRunVariance <- function(X,q){
  gamma0 <- autocov(X,0)
  LRV <- gamma0
  weights <- 1 - (1:q)/(q+1)
  for(i in 1:q){
    LRV <- LRV + 2 * weights[i] * autocov(X,i)
  }
  return(LRV)
}


compute_Gaussian_fast <- function(model,M){
  # This function computes mu_h, Phi_h, and Sigma_h that are s.t.
  # x_{t+h}+...+x_{t+1}|N(mu_h + Phi_h·x_t, Sigma_h),
  # where x_t follows the Gaussian VAR(1):
  # x_{t+1} = mu + Phi·x_t + Sigma^{1/2}·epsilon_{t+1},
  #          where epsilon_{t+1}~N(O,Id).
  # M determines (and decomposes) the maturities h of interest;
  # it is of dimension kx2. Let us denote its entries by m_{i,1} and
  # m_{i,2} (for row i), and the i'th horizon by h_{i}. The maturities are
  # obtained recursively using:
  #      h_{i+1} = h_{i-1}*m_{i,2} + h_{i-2}*m_{i,1}, with h_1=1 and h_0=0.
  # By convention, the first row of M is [1,0]. If, for instance:
  #      M[2,] = [4, 0] => h_2 = 4*1  + 0*0 = 4
  #      M[3,] = [3, 1] => h_3 = 3*4  + 1*1 = 13
  #      M[4,] = [4, 0] => h_4 = 4*13 + 0*4 = 52
  # This decomposition speeds up the calculation.

  mu  <- model$mu
  Phi <- model$Phi
  Sigma <- model$Sigma # covariance matrix

  n <- dim(Phi)[1]

  Id_n  <- diag(n)
  Id_n2 <- diag(n*n)
  Id_Phi_1 <- solve(Id_n - Phi)

  PhiPhi <- Phi %x% Phi
  Id_PhiPhi_1 <- solve(Id_n2 - PhiPhi)

  Sigma_star <- Id_Phi_1 %*% Sigma %*% t(Id_Phi_1)

  # Compute first column of A:
  Phi_exp_h_1 <- Phi
  Phi_exp_h_2 <- diag(n)

  PhiPhi_exp_h_1 <- PhiPhi
  PhiPhi_exp_h_2 <- diag(n*n)

  h_1 <- 1
  h_2 <- 0

  MU    <- array(NaN,c(n,1,dim(M)[1]))
  PHI   <- array(NaN,c(n,n,dim(M)[1]))
  SIGMA <- array(NaN,c(n,n,dim(M)[1]))
  H <- NULL

  for(i in 1:dim(M)[1]){
    h <- h_1*M[i,1] + h_2*M[i,2]

    Phi_exp_h <- (Phi_exp_h_1 %^% M[i,1]) %*% (Phi_exp_h_2 %^% M[i,2])
    PhiPhi_exp_h <- (PhiPhi_exp_h_1 %^% M[i,1]) %*% (PhiPhi_exp_h_2 %^% M[i,2])

    K_h   <- Id_Phi_1 %*% (Id_n - Phi_exp_h)
    mu_h  <- Id_Phi_1 %*% (h*Id_n - Phi %*% K_h) %*% mu
    Phi_h <- Phi %*% K_h
    aux <- Phi %*% K_h %*% Sigma_star

    Aux <- PhiPhi %*% Id_PhiPhi_1 %*% (Id_n2 - PhiPhi_exp_h) %*% c(Sigma_star)
    vecSigma_h <- h*c(Sigma_star) - c(aux) - c(t(aux)) + Aux

    # For next iteration:
    Phi_exp_h_2 <- Phi_exp_h_1
    Phi_exp_h_1 <- Phi_exp_h
    PhiPhi_exp_h_2 <- PhiPhi_exp_h_1
    PhiPhi_exp_h_1 <- PhiPhi_exp_h
    h_2 <- h_1
    h_1 <- h

    MU[,,i]  <- mu_h
    PHI[,,i] <- Phi_h
    SIGMA[,,i] <- matrix(vecSigma_h,n,n)
    H <- c(H,h)
  }
  return(list(MU=MU,PHI=PHI,SIGMA=SIGMA,H=H))
}

