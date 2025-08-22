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
#
# par(mfrow=c(2,1))
# plot(z[,1],type="l",lwd=2)
# lines(res_KH$ksi_matrix[,1],col="red")
# plot(z[,2],type="l",lwd=2)
# lines(res_KH$ksi_matrix[,2],col="red")
