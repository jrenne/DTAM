
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
  invisible(NULL)
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


compute_Gaussian_fast <- function(model,Maturities_decompo){
  # This function computes mu_h, Phi_h, and Sigma_h that are s.t.
  # x_{t+h}+...+x_{t+1}|N(mu_h + Phi_h·x_t, Sigma_h),
  # where x_t follows the Gaussian VAR(1):
  # x_{t+1} = mu + Phi·x_t + Sigma^{1/2}·epsilon_{t+1},
  #          where epsilon_{t+1}~N(O,Id).
  # Maturities_decompo determines (and decomposes)
  # these maturities; it is of dimension k x k. The k-th row determines
  # how the k-th horizon will be decomposed into the previous ones. Denote
  # the (i,j) entry of Maturities_decompo with m(i,j). Denote the i-th horizon
  # with h(i). We have:
  # h(i) = m(i,1)*1 + m(i,2)*h(1) + ... + m(i,i-1)*h(i-1).
  # For instance:
  #   Maturities_decompo[1,] = [2, 0, 0] => h_1 = 2
  #   Maturities_decompo[2,] = [1, 2, 0] => h_2 = 1*1 + 2*2 = 5
  #   Maturities_decompo[3,] = [2, 0, 2] => h_3 = 2*1 + 0*2 + 2*5 = 12.
  # This decomposition speeds up the calculation by decreasing the number of
  # matrix multiplications.

  mu  <- model$mu
  Phi <- model$Phi
  Sigma <- model$Sigma # covariance matrix

  n <- dim(Phi)[1]

  k <- dim(Maturities_decompo)[1]

  Id_n  <- diag(n)
  Id_n2 <- diag(n*n)
  Id_Phi_1 <- solve(Id_n - Phi)

  PhiPhi <- Phi %x% Phi
  Id_PhiPhi_1 <- solve(Id_n2 - PhiPhi)

  Sigma_star <- Id_Phi_1 %*% Sigma %*% t(Id_Phi_1)

  MU    <- array(NaN,c(n,1,k))
  PHI   <- array(NaN,c(n,n,k))
  SIGMA <- array(NaN,c(n,n,k))
  H <- NULL

  all_Phi_h <- array(NaN,c(n,n,k))
  all_PhiPhi_h <- array(NaN,c(n*n,n*n,k))

  for(i in 1:k){
    Phi_exp_h    <- Id_n
    PhiPhi_exp_h <- Id_n2
    h_i <- 0
    for(j in 1:i){
      if(j==1){
        Aux <- Phi
        AAux <- PhiPhi
        aux <- 1
      }else{
        Aux <- all_Phi_h[,,j-1]
        AAux <- all_PhiPhi_h[,,j-1]
        aux <- H[j-1]
      }
      Phi_exp_h    <- Phi_exp_h    %*% (Aux %^% Maturities_decompo[i,j])
      PhiPhi_exp_h <- PhiPhi_exp_h %*% (AAux %^% Maturities_decompo[i,j])
      h_i  <- h_i + aux*Maturities_decompo[i,j]
    }
    all_Phi_h[,,i]    <- Phi_exp_h
    all_PhiPhi_h[,,i] <- PhiPhi_exp_h

    K_h   <- Id_Phi_1 %*% (Id_n - Phi_exp_h)
    mu_h  <- Id_Phi_1 %*% (h_i*Id_n - Phi %*% K_h) %*% mu
    Phi_h <- Phi %*% K_h
    aux <- Phi %*% K_h %*% Sigma_star

    Aux <- PhiPhi %*% Id_PhiPhi_1 %*% (Id_n2 - PhiPhi_exp_h) %*% c(Sigma_star)
    vecSigma_h <- h_i*c(Sigma_star) - c(aux) - c(t(aux)) + Aux

    # Store results:
    MU[,,i]  <- mu_h
    PHI[,,i] <- Phi_h
    SIGMA[,,i] <- matrix(vecSigma_h,n,n)
    H <- c(H,h_i)
  }
  return(list(MU=MU,PHI=PHI,SIGMA=SIGMA,H=H))
}

# Maturities_decompo <- diag(4)
# Maturities_decompo[2,1:2] <- c(0,4)
# Maturities_decompo[3,1:3] <- c(0,1,3)
# Maturities_decompo[4,1:4] <- c(1,1,1,4)
#
# mu <- matrix(c(0,1),2,1)
# Phi <- .9*diag(2)
# Phi[2,1] <- .4
# Sigma <- .1*diag(2)
# Sigma[1,2] <- .005
# Sigma[2,1] <- .005
#
# model <- list(mu=mu,Phi=Phi,Sigma=Sigma)
#
# RES <- compute_Gaussian_fast(model,Maturities_decompo)
# maxH <- RES$H[length(RES$H)]
#
# u <- matrix(1,2,1)
# res <- reverse.MHLT(psi.GaussianVAR,u1 = u,H = maxH,psi.parameterization = model)
#
# m <- dim(Maturities_decompo)[1]
# mu_h <- matrix(RES$MU[,,m],ncol=1)
# Phi_h <- RES$PHI[,,m]
# Sigma_h <- RES$SIGMA[,,m]
#
# t(u) %*% mu_h + (t(u) %*% Sigma_h %*% u)/2
# t(u) %*% Phi_h

simul.GVAR_OLD <- function(model,nb.sim,x0=NaN){
  n <- dim(model$Phi)[1]
  if(is.na(x0[1])){
    x0 <- solve(diag(n) - model$Phi) %*% model$mu
  }

  if(is.null(model$Sigma12)){
    Sigma12 <- t(chol(model$Sigma))
  }else{
    Sigma12 <- model$Sigma12
  }

  X <- c(x0)
  x <- x0
  for(t in 2:nb.sim){
    x <- model$mu + model$Phi %*% x + Sigma12 %*% rnorm(dim(Sigma12)[2])
    X <- rbind(X,c(x))
  }
  return(X)
}

simul.GVAR <- function(model,nb.sim,x0=NaN,nb.replic=1){
  # Simulate a Gaussian VAR model.
  # If nb.replic > 1, simluate several processes of length nb.sim in parallel;
  #.    the output matrix then is an array (three dimensions).

  n <- dim(model$Phi)[1]
  if(is.na(x0[1])){
    x0 <- solve(diag(n) - model$Phi) %*% model$mu
  }

  if(is.null(model$Sigma12)){
    Sigma12 <- t(chol(model$Sigma))
  }else{
    Sigma12 <- model$Sigma12
  }

  if(nb.replic==1){
    X <- c(x0)
    x <- x0
    for(t in 2:nb.sim){
      x <- model$mu + model$Phi %*% x + Sigma12 %*% rnorm(dim(Sigma12)[2])
      X <- rbind(X,c(x))
    }
  }else{
    # In that case, we simulate nb.replic paths in parallel
    X <- array(NaN,c(nb.sim,n,nb.replic))
    X[1,,] <- x0
    x <- matrix(x0,n,nb.replic)
    for(t in 2:nb.sim){
      x <- matrix(model$mu,n,nb.replic) + model$Phi %*% x + Sigma12 %*% matrix(rnorm(dim(Sigma12)[2]*nb.replic),dim(Sigma12)[2],nb.replic)
      X[t,,] <- x
    }
  }
  return(X)
}

plot_isodensity_gaussian <- function(Ew, Vw, n = 200,
                                     min_w = NA_real_, max_w = NA_real_,
                                     prob_levels = c(0.50, 0.80, 0.95, 0.99),
                                     n_std = 4,
                                     maintitle = "") {
  stopifnot(length(Ew) == 2, all(dim(Vw) == c(2, 2)))

  sds <- sqrt(diag(Vw))

  # Default axis range (if not supplied)
  if (!is.finite(min_w)) min_w <- Ew[1] - n_std * sds[1]
  if (!is.finite(max_w)) max_w <- Ew[1] + n_std * sds[1]

  # Grid for density
  x1 <- seq(Ew[1] - n_std * sds[1], Ew[1] + n_std * sds[1], length.out = n)
  x2 <- seq(Ew[2] - n_std * sds[2], Ew[2] + n_std * sds[2], length.out = n)

  X1 <- matrix(rep(x1, each = n), nrow = n)
  X2 <- matrix(rep(x2, times = n), nrow = n)

  detV <- det(Vw)
  invV <- solve(Vw)

  dX1 <- X1 - Ew[1]
  dX2 <- X2 - Ew[2]

  # Quadratic form Q = (w-Ew)' V^{-1} (w-Ew)
  Q <- invV[1,1]*dX1^2 + 2*invV[1,2]*dX1*dX2 + invV[2,2]*dX2^2

  # Density
  logdens <- -0.5 * (2 * log(2*pi) + log(detV) + Q)
  dens <- exp(logdens)

  # Convert probability-mass levels to density levels
  # In 2D Gaussian: Q ~ chi^2(df=2) on ellipses
  Q_levels <- qchisq(prob_levels, df = 2)
  dens_levels <- exp(-0.5 * Q_levels) / (2*pi*sqrt(detV))

  contour(x1, x2, dens,
          las = 1,
          xlim = c(min_w, max_w),
          ylim = c(min_w, max_w),
          levels = dens_levels,
          labels = paste0(100 * prob_levels, "%"),
          drawlabels = TRUE,
          xlab = expression(r[t]),
          ylab = expression(x[t]),
          main = maintitle,
          labcex = 0.8,
          vfont = c("sans serif", "bold"),
          lwd = 2, col = "grey")

  grid()

  # Mark unconditional mean
  points(Ew[1], Ew[2], pch = 3, lwd = 2)
  text(Ew[1], Ew[2], labels = expression(E[t](w[T])), pos = 4)

  invisible(list(x1 = x1, x2 = x2, dens = dens,
                 dens_levels = dens_levels, prob_levels = prob_levels))
}


