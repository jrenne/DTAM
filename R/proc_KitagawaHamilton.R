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
  #  If it is a 3-dimensional array, then it means that the transition
  #  probabilities vary (deterministically) over time. The third dimension
  #  has to be T.
  # Eta is a matrix of dimension T x J, its (t,j) entry is the evaluation of
  #  the pdf of the distribution F_{t}|z_{t}=e_j,I_{t-1}.
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
  if(length(dim(Omega))==2){
    Omega_2h <- t(Omega)
  }else{
    Omega_2h <- t(Omega[,,1])
  }
  for (i in 1:10){
    Omega_2h <- Omega_2h %*% Omega_2h
  }
  stat_distri[] <- Omega_2h[, 1]
  ksi_t[] <- stat_distri

  loglik_vec <- NULL

  for (t in 1:nb_dates){

    if(length(dim(Omega))==2){
      Omega_t <- Omega
    }else{# In that case, Omega changes over time
      Omega_t <- Omega[,,t]
    }
    eta_t <- matrix(Eta[t, ], ncol = 1)
    # When all entries of eta_t are small and negative, the risk is that
    # exp(eta_t) is 0 for some entries. To avoid numerical instability,
    # we first shift all values of eta_t such that that the max one is
    # normalized at 0. "cst" is further reintroduced when needed.
    cst <- max(eta_t)
    eta_t_minus_max <- eta_t - cst
    exp_eta_t_minus_max <- exp(eta_t_minus_max)
    ksi_t <- (t(Omega_t) %*% ksi_t) * exp_eta_t_minus_max
    ksi_t[ksi_t==0] <- 10^(-60)
    loglik_vec <- c(loglik_vec, log(sum(ksi_t)) + cst)
    normalisation_factor <- sum(ksi_t)
    ksi_t <- ksi_t/normalisation_factor
    ksi_matrix[t, ] <- t(ksi_t)
  }

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

    if(length(dim(Omega))==2){
      Omega_t <- Omega
    }else{# In that case, Omega changes over time
      Omega_t <- Omega[,,t]
    }

    ksi_tt <- matrix(res_filter$ksi_matrix[t,],ncol=1)
    AUX <- (((ksi_tt %*% t(vec1)) * Omega_t)/((vec1 %*%
                                                 t(ksi_tt)) %*% Omega_t))
    AUX[is.na(AUX)] <- 0
    ksi_tT[t, ] <- AUX %*% matrix(ksi_tT[t + 1,], ncol = 1)
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
  # Pi is the matrix of transition probabilities (with rows that sum to one)
  # This procedure constructs the matrix  of transition probabilities of
  # z_{t} %x% z_{t-1}.
  # Pi can be an array of dimension J x J x T if the matrix of transition
  # probabilities depends on time.

  J <- dim(Pi)[1]
  Id <- diag(J)

  if(length(dim(Pi))==2){# Pi is a matrix in this case
    PI <- NULL
    for(i in 1:J){
      e_i <- matrix(Id[i,],nrow=1)
      aux <- matrix(Pi[i,],nrow=1) %x% e_i
      PI <- rbind(PI,aux)
    }
    PI <- PI %x% matrix(1,J,1)
  }else{
    # Pi is an array (dimension 3)
    PI <- array(0,c(J^2,J^2,dim(Pi)[3]))
    for(i in 1:J){
      for(j in 1:J){
        PI[((i-1)*J+1):(i*J),(j-1)*J+i,] <- Pi[i,j,] %x% matrix(1,J,1)
      }
    }
  }
  return(PI)
}



compute_LT_RS <- function(alpha,Pi,Maturities_decompo,
                          indic_add_current=FALSE){
  # This function computes A_h s.t.
  # E_t[exp(alpha'(z_{t+1}+...+z_{t+h}))] = [1,...,1]·A_h'·z_t,
  # where z_t follows a J-regime Markov-Switching process with matrix of
  # transition probabilities Pi (rows sum to one).
  # The A_h vector are collected in matrix A (J*k), where k is the number of
  # maturities of interest. Maturities_decompo determines (and decomposes)
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
  # ----------------------------------------------------------------------------
  # If "indic_add_current" is TRUE, then what is actually computed is
  #     exp(alpha'(z_{t}) * E_t[exp(alpha'(z_{t+1}+...+z_{t+h}))]
  # ----------------------------------------------------------------------------

  J <- dim(Pi)[1]
  k <- dim(Maturities_decompo)[1]

  vec1 <- matrix(1,J,1)

  all_DP <- array(NaN,c(J,J,k))

  A <- NULL
  H <- NULL
  DP1 <- diag(c(exp(alpha))) %*% t(Pi)
  for(i in 1:k){
    DP_i <- diag(J)
    h_i <- 0
    for(j in 1:i){
      if(j==1){
        Aux <- DP1
        aux <- 1
      }else{
        Aux <- all_DP[,,j-1]
        aux <- H[j-1]
      }
      DP_i <- DP_i %*% (Aux %^% Maturities_decompo[i,j])
      h_i  <- h_i + aux*Maturities_decompo[i,j]
    }
    all_DP[,,i] <- DP_i
    A_i <- t(DP_i) %*% vec1
    A <- cbind(A,A_i)
    H <- c(H,h_i)
  }

  if(indic_add_current){
    A <- diag(c(exp(alpha))) %*% A
  }

  return(list(A=A,H=H))
}


# library(mvtnorm)
# library(expm)
# library(DTAM)
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
# Omega <- array(Omega,c(3,3,TT))
#
# res_KH     <- KH_filter(Omega, Eta = f_Eta(F,M,N))
# #res_KH_new <- KH_filter_new(Omega, Eta = f_Eta(F,M,N))
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
#
# alpha <- c(-.03,-.01,-.05)
# Maturities_decompo <- diag(4)
# Maturities_decompo[2,1:2] <- c(0,4)
# Maturities_decompo[3,1:3] <- c(0,1,3)
# Maturities_decompo[4,1:4] <- c(1,1,1,4)
#
# Pi <- matrix(Omega,J,J)
# res <- compute_LT_RS(alpha,Pi,Maturities_decompo,indic_add_current = FALSE)
#
#
# D <- diag(exp(alpha))
#
# print(res$A)
# print(matrix(1,1,3) %*% ((D %*% t(Omega))%^%max(res$H)))
#
#
#
