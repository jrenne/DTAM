KH_filter <- function(Omega, Eta){

  # Model: ---------------------------------------------------------------------
  # F_{t}|z_{t},I_{t-1} ~ f(z_{t},F_{t-1}),
  #  where I_{t-1} is the information available on date t-1.
  # F_{t} is of dimension n x 1.
  # ----------------------------------------------------------------------------
  # z_t follows a Markov chain with matrix of transition probabilities Omega,
  # Omega (of dim. J x J) is such that the entries of each of its rows sum to 1.
  # ----------------------------------------------------------------------------
  # F is a matrix of dimension T x n.
  # Omega is a matrix of dimension J x J (rows sum to 1).
  # Eta is a matrix of dimension T x J, its (t,j) entry is the evaluation of
  #       the pdf of the distribution F_{t}|z_{t}=e_j,I_{t-1}.
  # ----------------------------------------------------------------------------

  J         <- dim(Omega)[1]  # number of regimes
  nb_dates  <- nrow(Eta)

  ksi_matrix   <- matrix(0, nrow = nb_dates, ncol = J)
  ksi_1_matrix <- matrix(0, nrow = nb_dates, ncol = J)
  ksi_t        <- matrix(0, nrow = J, ncol = 1)

  vec_1_J     <- matrix(1, nrow = J, ncol = 1)

  # Compute stationary distribution:
  Omega_2h <- matrix(0, nrow = J, ncol = J)
  stat_distri <- matrix(0, nrow = J, ncol = 1)
  Omega_2h <- t(Omega)
  for (i in 1:10){
    Omega_2h <- Omega_2h %*% Omega_2h
  }
  stat_distri[] <- Omega_2h[, 1]
  ksi_t[] <- stat_distri

  for (t in 1:nb_dates){
    ksi_1_matrix[t, ] <- t(ksi_t)  # will be used to compute log-likelihood

    eta_t <- matrix(Eta[t, ],ncol=1)  # eta_matrix is LOG-likelihood
    eta_t <- exp(eta_t)
    # Use update formula:
    ksi_t <- (t(Omega) %*% ksi_t) * eta_t
    normalisation_factor <- sum(ksi_t)
    ksi_t <- ksi_t / normalisation_factor

    ksi_matrix[t, ] <- t(ksi_t)
  }

  # Compute log-likelihood:
  loglik_mat <- (ksi_1_matrix %*% Omega) * exp(Eta)
  loglik_vec <- log(loglik_mat %*% vec_1_J)
  loglik <- sum(loglik_vec)

  return(list(ksi_matrix = ksi_matrix,
              loglik_vec = loglik_vec,
              loglik = loglik))
}

KH_smoother <- function(Omega, Eta){
  res_filter <- KH_filter(Omega, Eta)

  nb_dates <- dim(Eta)[1]
  J <- dim(Eta)[2]
  vec1 <- matrix(1,J,1)

  ksi_tT <- res_filter$ksi_matrix

  for(t in (nb_dates-1):1){
    ksi_tt <- matrix(res_filter$ksi_matrix[t,],ncol=1)
    ksi_tT[t,] <- (((ksi_tt %*% t(vec1))*Omega)/((vec1 %*% t(ksi_tt))%*%Omega)) %*%
      matrix(ksi_tT[t+1,],ncol=1)
  }
  return(ksi_tT)
}

f_Eta <- function(F,M,N){
  # Compute log-likelihood conditional on each regime:

  J         <- ncol(M)  # number of regimes
  nb_dates  <- nrow(F)
  nb_obsvar <- ncol(F)

  vec_1_dates <- matrix(1, nrow = nb_dates, ncol = 1)
  eta_matrix  <- matrix(0, nrow = nb_dates, ncol = J)

  variances4regime <- matrix(0, nrow = nb_obsvar, ncol = 1)
  Covariance       <- matrix(0, nrow = nb_obsvar, ncol = nb_obsvar)
  for (j in 1:J){
    variances4regime <- (N[, j])^2
    Covariance <- diag(variances4regime)
    eta_matrix[, j] <- log(dmvnorm(F - vec_1_dates %*% t(M[, j]), sigma=Covariance))
  }
  return(eta_matrix)
}



simul_RS <- function(Omega,TT,ini_state = NaN){
  # Omega is the matrix of transition probabilities, its rows sum to 1.

  # Number of regimes:
  J = dim(Omega)[1]

  if(is.na(ini_state[1])){
    # in this case, the starting state is the most likely (unconditionally):
    # For that, compute stationary distribution:
    Omega_2h <- matrix(0, nrow = J, ncol = J)
    stat_distri <- matrix(0, nrow = J, ncol = 1)
    Omega_2h <- t(Omega)
    for (i in 1:10){
      Omega_2h <- Omega_2h %*% Omega_2h
    }
    stat_distri[] <- Omega_2h[, 1]
    ini_state <- which(stat_distri==max(stat_distri))[1]
  }

  state <- ini_state

  states <- matrix(0,TT,J)
  states[1,state] <- 1
  for(t in 2:TT){
    p    <- Omega[state,]
    cump <- cumsum(p)
    u    <- runif(1)
    state <- which(cump > u)[1]
    states[t,state] <- 1
  }

  return(states)
}


make_Pi_z_kron_z <- function(Pi){
  # Pi is the matrix of transition probabilties (with rows that sum to one)
  # This procedure constructs the matrix  of transition probabilities of
  # z_{t} %x% z_{t-1}.
  J <- dim(Pi)[1]
  Id <- diag(J)

  PI <- NULL
  for(i in 1:J){
    e_i <- matrix(Id[i,],nrow=1)
    aux <- matrix(Pi[i,],nrow=1) %x% e_i
    PI <- rbind(PI,aux)
  }
  PI <- PI %x% matrix(1,J,1)
  return(PI)
}



compute_LT_RS <- function(alpha,Pi,M){
  # This function computes A_h s.t.
  # E_t[exp(alpha'(z_{t+1}+...+z_{t+h}))] = [1,...,1]·A_h'·z_t,
  # where z_t follows a J-regime Markov-Switching process with matrix of
  # transition probabilities Pi (rows sum to one).
  # The A_h vector are collected in matrix A (J*k), where k is the number of
  # maturities of interest. M determines (and decomposes) these maturities;
  # it is of dimension kx2. Let us denote its entries by m_{i,1} and
  # m_{i,2} (for row i), and the i'th horizon by h_{i}. The maturities are
  # obtained recursively using:
  #      h_{i+1} = h_{i-1}*m_{i,2} + h_{i-2}*m_{i,1}, with h_1=1 and h_0=0.
  # By convention, the first row of M is [1,0]. If, for instance:
  #      M[2,] = [4, 0] => h_2 = 4*1  + 0*0 = 4
  #      M[3,] = [3, 1] => h_3 = 3*4  + 1*1 = 13
  #      M[4,] = [4, 0] => h_4 = 4*13 + 0*4 = 52
  # This decomposition speeds up the calculation.

  J <- dim(Pi)[1]
  vec1 <- matrix(1,J,1)

  # Compute first column of A:
  DP_h_1 <- diag(c(exp(alpha))) %*% t(Pi)
  DP_h_2 <- diag(J)
  A <- NULL
  for(i in 1:dim(M)[1]){
    DP_h <- (DP_h_1 %^% M[i,1]) %*% (DP_h_2 %^% M[i,2])
    DP_h_2 <- DP_h_1 # for next iteration
    DP_h_1 <- DP_h   # for next iteration
    A_h <- t(DP_h) %*% vec1
    A <- cbind(A,A_h)
  }
  return(A)
}

#
# Omega <- matrix(c(.8,.1,0,.2,.8,.2,0,.1,.8),3,3)
# J <- dim(Omega)[1]
# TT <- 200
# z <- simul_RS(Omega,TT=TT,ini_state = NaN)
#
# nF <- 2
# M <- matrix(c(0,1,1,0,1,1),nF,J)
# N <- .5*matrix(c(.5,1,.3,1,2,.6),nF,J)
#
# F <- z %*% t(M) + (z %*% t(N))*matrix(rnorm(nF*TT),TT,nF)
#
# res_KH <- KH_filter(Omega, Eta = f_Eta(F,M,N))
#
# res_smoother <- res_smoother <- KH_smoother(Omega, Eta = f_Eta(F,M,N))
#
# par(mfrow=c(2,1))
# plot(z[,1],type="l",lwd=2)
# lines(res_KH$ksi_matrix[,1],col="red")
# lines(res_smoother[,1],col="blue")
# plot(z[,2],type="l",lwd=2)
# lines(res_KH$ksi_matrix[,2],col="red")
# lines(res_smoother[,2],col="blue")
#
# alpha <- c(-.03,-.01,-.05)
# M <- t(matrix(c(1,0,4,0,3,1,4,0),2,4))
#
# A <- compute_LT_RS(alpha,Omega,M)
#
# D <- diag(exp(alpha))
#
# print(A)
# print(matrix(1,1,3) %*% ((D %*% t(Omega))%^%52))
#
#
# Phi <- matrix(1:4,2,2)
#
#
# solve(diag(2) - Phi) %*% (diag(2) - Phi%^%3)
# (diag(2) - Phi%^%3) %*% solve(diag(2) - Phi)
#
# t(solve(diag(2) - Phi))
# solve(diag(2) - t(Phi))
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
# RES <- compute_Gaussian_fast(model,M)
# maxH <- RES$H[length(RES$H)]
#
# u <- matrix(1,2,1)
# res <- reverse.MHLT(psi.GaussianVAR,u1 = u,H = maxH,psi.parameterization = model)
#
# m <- dim(M)[1]
# mu_h <- matrix(RES$MU[,,m],ncol=1)
# Phi_h <- RES$PHI[,,m]
# Sigma_h <- RES$SIGMA[,,m]
#
# t(u) %*% mu_h + (t(u) %*% Sigma_h %*% u)/2
# t(u) %*% Phi_h
#
#
