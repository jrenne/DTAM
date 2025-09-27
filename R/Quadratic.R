

vecX <- function(X) {
  p <- nrow(X)
  q <- ncol(X)
  vecX <- matrix(0, p * q, 1)

  for (j in 1:q) {
    for (i in 1:p) {
      vecX[j * p - p + i, 1] <- X[i, j]
    }
  }
  return(vecX)
}

vec_1 <- function(vecX) {
  n2 <- nrow(vecX)
  n <- sqrt(n2)
  X <- matrix(0, n, n)

  for (j in 1:n) {
    for (i in 1:n) {
      X[i, j] <- vecX[j * n - n + i, 1]
    }
  }
  return(X)
}

make_permut_matrices <- function(n) {
  n2 <- n * n
  ntot <- n + n2
  ntot2 <- ntot * ntot

  I_A <- matrix(0, ntot2, n2)
  I_B <- matrix(0, ntot2, n2 * n)
  I_C <- matrix(0, ntot2, n * n2)
  I_D <- matrix(0, ntot2, n2 * n2)

  AUXA <- matrix(1, n, n)
  AUXB <- matrix(1, n2, n)
  AUXC <- matrix(1, n, n2)
  AUXD <- matrix(1, n2, n2)

  AUXAUXA <- matrix(0, ntot, ntot)
  AUXAUXB <- matrix(0, ntot, ntot)
  AUXAUXC <- matrix(0, ntot, ntot)
  AUXAUXD <- matrix(0, ntot, ntot)

  AUXAUXA[1:n, 1:n] <- AUXA
  AUXAUXB[(n + 1):ntot, 1:n] <- AUXB
  AUXAUXC[1:n, (n + 1):ntot] <- AUXC
  AUXAUXD[(n + 1):ntot, (n + 1):ntot] <- AUXD

  I <- diag(ntot2)
  count <- 1
  count_A <- 1
  count_B <- 1
  count_C <- 1
  count_D <- 1

  for (i in 1:ntot) {
    for (j in 1:ntot) {
      if (AUXAUXA[i, j] == 1) {
        I_A[, count_A] <- I[, count]
        count_A <- count_A + 1
      }
      if (AUXAUXB[i, j] == 1) {
        I_B[, count_B] <- I[, count]
        count_B <- count_B + 1
      }
      if (AUXAUXC[i, j] == 1) {
        I_C[, count_C] <- I[, count]
        count_C <- count_C + 1
      }
      if (AUXAUXD[i, j] == 1) {
        I_D[, count_D] <- I[, count]
        count_D <- count_D + 1
      }
      count <- count + 1
    }
  }

  return(list(I_A = I_A, I_B = I_B, I_C = I_C, I_D = I_D))
}


vcov_matrix_Quadratic <- function(mu, Phi, Sigma, BigMats=make_BigMats(dim(Phi)[1])) {
  # The model is:
  # x_t = mu + Phi·x_{t-1} + eps_t, with eps_t ~ N(0,Sigma).
  # The extended state vector is:
  # w_t = [x_t',vec(x_t·x_t')]' or
  # w_t = [x_t',vech(x_t·x_t')]' (the outputs then have the "_red" suffix)

  n <- nrow(Phi)
  n2 <- n * n
  n4 <- n2 * n2
  nnp12 <- n * (n + 1) / 2
  n_total <- n + nnp12
  n_total2 <- n_total * n_total

  M_vec_vech <- BigMats$M_vec_vech
  M_vech_vec <- BigMats$M_vech_vec

  Id <- diag(n_total)
  Idn <- diag(n)
  Idn2 <- diag(n2)
  Idnn2 <- diag(n + n2)

  M_YZ <- matrix(0, n_total, n + n2)
  M_YZ[1:n, 1:n] <- Idn
  M_YZ[(n + 1):n_total, (n + 1):(n + n2)] <- M_vech_vec

  M_ZY <- matrix(0, n + n2, n_total)
  M_ZY[1:n, 1:n] <- Idn
  M_ZY[(n + 1):(n + n2), (n + 1):n_total] <- M_vec_vech

  Phi_tilde <- matrix(0, n + n2, n + n2)
  Phi_tilde[1:n, 1:n] <- Phi
  Phi_tilde[(n + 1):(n + n2), 1:n] <-
    (mu %x% Phi) + (Phi %x% mu)
  Phi_tilde[(n + 1):(n + n2), (n + 1):(n + n2)] <- (Phi %x% Phi)

  Phi_tilde_red <- M_YZ %*% Phi_tilde %*% M_ZY

  mu_tilde <- matrix(0, n + n2, 1)
  mu_tilde[1:n, ] <- mu
  aux <- mu %*% t(mu) + Sigma
  mu_tilde[(n + 1):(n + n2), ] <- vecX(aux)

  mu_tilde_red <- M_YZ %*% mu_tilde

  E <- solve(diag(n + n2) - Phi_tilde) %*% mu_tilde

  Ex <- solve(diag(n) - Phi) %*% mu

  Gamma <- matrix(0, n2, n)
  cond_exp <- mu + Phi %*% Ex
  Gamma <- (diag(n) %x% cond_exp) + (cond_exp %x% diag(n))

  Omega_tilde <- matrix(0, n_total, n_total)
  Omega_tilde[1:n, 1:n] <- Sigma
  Omega_tilde[(n + 1):n_total, 1:n] <-
    M_vech_vec %*% Gamma %*% Sigma
  Omega_tilde[1:n, (n + 1):n_total] <-
    Sigma %*% t(Gamma) %*% t(M_vech_vec)

  Aux1 <- BigMats$Aux1
  Aux2 <- BigMats$Aux2

  aux3 <- (vecX(Sigma) %x% diag(n2))
  Phi2_tilde <- Phi_tilde[(n + 1):(n + n2), ]
  mumu <- Kronecker_matmat(mu, mu)
  aux4 <- mumu + Phi2_tilde %*% E
  aux6 <- vecX(Sigma %x% Sigma)

  vecSigma22 <- Aux1 %*% aux3 %*% aux4 + Aux2 %*% aux6
  Sigma22 <- matrix(0, nnp12, nnp12)

  for (i in 1:nnp12) {
    Sigma22[, i] <- vecSigma22[(i - 1) * nnp12 + (1:nnp12)]
  }

  Omega_tilde[(n + 1):n_total, (n + 1):n_total] <- Sigma22

  V_red <- Omega_tilde
  for (i in 1:100) {
    V_red <- Omega_tilde + Phi_tilde_red %*% V_red %*% t(Phi_tilde_red)
  }

  V <- M_ZY %*% V_red %*% t(M_ZY)

  return(list(
    mu_tilde = mu_tilde,
    Phi_tilde = Phi_tilde,
    V_red = V_red,
    E = E,
    V = V
  ))
}

make_BigMats <- function(n) {
  n2 <- n * n
  n4 <- n2 * n2
  nM <- n * (n + 1) / 2
  nM2 <- nM * nM

  Idn <- diag(n)
  Idn2 <- diag(n2)

  M_vech_vec <- matrix(0, nM, n2)
  M_vec_vech <- matrix(0, n2, nM)

  Lambda_n <- matrix(0, n2, n2)
  for (i in 1:n) {
    e_i <- Idn[, i]
    for (j in 1:n) {
      e_j <- Idn[, j]
      Lambda_n[(i - 1) * n + (1:n), (j - 1) * n + (1:n)] <- e_j %*% t(e_i)
    }
  }

  Aux <- diag(n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      Aux[i, j] <- 1
    }
  }

  vec_Aux <- vecX(Aux)
  count <- 1
  for (i in 1:n2) {
    if (vec_Aux[i, 1] == 1) {
      M_vech_vec[count, ] <- Idn2[i, ]
      count <- count + 1
    }
  }

  Lambda_n_minus_diag <- Lambda_n
  diag(Lambda_n_minus_diag) <- 0
  M_vec_vech <- (Idn2 + Lambda_n_minus_diag) %*% t(M_vech_vec)

  MkronM <- M_vech_vec %x% M_vech_vec

  aux1 <- (Idn2 + Lambda_n) %x% (Idn2 + Lambda_n)
  aux2 <- diag(n) %x% (Lambda_n %x% diag(n))
  Aux1 <- MkronM %*% aux1 %*% aux2

  Aux2 <- MkronM %*% (Idn2 %x% (Idn2 + Lambda_n))

  return(list(
    Aux1 = Aux1,
    Aux2 = Aux2,
    Lambda_n = Lambda_n,
    M_vech_vec = M_vech_vec,
    M_vec_vech = M_vec_vech
  ))
}

make_matrices_cond_mean_variance_Quadratic <- function(model, indic_compute_V) {
  # The model is:
  # x_t = mu + Phi·x_{t-1} + eps_t, with eps_t ~ N(0,Sigma).
  # The extended state vector is:
  # w_t = [x_t',vec(x_t·x_t')]' or
  # w_t = [x_t',vech(x_t·x_t')]' (the outputs then have the "_red" suffix)
  # The function returns vectors and matrices (mu_tilde, Phi_tilde, nu, and Psi)
  #   that are such that:
  #     E_t(w_{t+1}) = mu_tilde + Phi_tilde·w_t
  #     Var_t(w_{t+1}) = nu + Psi·w_t

  mu <- model$mu
  Phi <- model$Phi
  Sigma <- model$Sigma

  n <- nrow(Phi)
  n2 <- n * n
  n4 <- n2 * n2
  ntot <- n + n2
  nnp12 <- n * (n + 1) / 2
  n_total <- n + nnp12
  n_total2 <- n_total * n_total
  n2np12 <- nnp12 * nnp12

  # Precompute matrices from helper functions
  BigMats <- make_BigMats(n)
  I_matrices <- make_permut_matrices(n)

  M_vec_vech <- BigMats$M_vec_vech
  M_vech_vec <- BigMats$M_vech_vec
  Lambda_n <- BigMats$Lambda_n
  Aux1 <- BigMats$Aux1
  Aux2 <- BigMats$Aux2

  I_A <- I_matrices$I_A
  I_B <- I_matrices$I_B
  I_C <- I_matrices$I_C
  I_D <- I_matrices$I_D

  # Identity matrices
  Id <- diag(n_total)
  Idn <- diag(n)
  Idn2 <- diag(n2)
  Idnn2 <- diag(n + n2)
  Idn2n4 <- diag(ntot * ntot)

  # Construct M_YZ and M_ZY
  M_YZ <- matrix(0, n_total, n + n2)
  M_YZ[1:n, 1:n] <- Idn
  M_YZ[(n + 1):n_total, (n + 1):(n + n2)] <- M_vech_vec

  M_ZY <- matrix(0, n + n2, n_total)
  M_ZY[1:n, 1:n] <- Idn
  M_ZY[(n + 1):(n + n2), (n + 1):n_total] <- M_vec_vech

  # Construct Phi_tilde
  Phi_tilde <- matrix(0, n + n2, n + n2)
  Phi_tilde[1:n, 1:n] <- Phi
  Phi_tilde[(n + 1):(n + n2), 1:n] <-
    (mu %x% Phi) + (Phi %x% mu)
  Phi_tilde[(n + 1):(n + n2), (n + 1):(n + n2)] <- (Phi %x% Phi)

  # Reduced Phi_tilde
  Phi_tilde_red <- M_YZ %*% Phi_tilde %*% M_ZY

  # Construct mu_tilde
  mu_tilde <- matrix(0, n + n2, 1)
  mu_tilde[1:n, ] <- mu
  aux <- mu %*% t(mu) + Sigma
  mu_tilde[(n + 1):(n + n2), ] <- vecX(aux)

  # Reduced mu_tilde
  mu_tilde_red <- M_YZ %*% mu_tilde

  # Extract Phi1_tilde and Phi2_tilde
  Phi1_tilde <- Phi_tilde[1:n, 1:(n + n2)]
  Phi2_tilde <- Phi_tilde[(n + 1):(n + n2), 1:(n + n2)]

  # Compute vecSigma11_const
  vecSigma11_const <- vecX(Sigma)

  # Compute vecSigma12 components
  aux_vecSigma12 <-
    (Sigma %x% (Idn2 + Lambda_n)) %*% (vecX(Idn) %x% Idn)
  vecSigma12_const <- aux_vecSigma12 %*% mu
  vecSigma12_Z <- aux_vecSigma12 %*% Phi1_tilde

  # Compute vecSigma21 components
  aux_vecSigma21 <-
    ((Idn2 + Lambda_n) %x% Sigma) %*% (Idn %x% Lambda_n) %*% (vecX(Idn) %x% Idn)
  vecSigma21_const <- aux_vecSigma21 %*% mu
  vecSigma21_Z <- aux_vecSigma21 %*% Phi1_tilde

  # Compute vecSigma22 components
  aux_vecSigma22 <-
    ((Idn2 + Lambda_n) %x% (Idn2 + Lambda_n)) %*%
    (Idn %x% (Lambda_n %x% Idn)) %*% (vecX(Sigma) %x% Idn2)
  vecSigma22_const <-
    aux_vecSigma22 %*% (mu %x% mu) +
    (Idn2 %x% (Idn2 + Lambda_n)) %*% vecX(Sigma %x% Sigma)
  vecSigma22_Z <- aux_vecSigma22 %*% Phi2_tilde

  # Compute nu and Psi
  nu <- I_A %*% vecSigma11_const +
    I_B %*% vecSigma21_const +
    I_C %*% vecSigma12_const +
    I_D %*% vecSigma22_const

  Psi <- I_B %*% vecSigma21_Z +
    I_C %*% vecSigma12_Z +
    I_D %*% vecSigma22_Z

  # Reduced nu and Psi
  nu_red <- (M_YZ %x% M_YZ) %*% nu
  Psi_red <- (M_YZ %x% M_YZ) %*% Psi %*% M_ZY

  # Compute E and E_red
  E <- solve(diag(n + n2) - Phi_tilde) %*% mu_tilde
  E_red <- M_YZ %*% E

  # Initialize V and V_red
  V <- matrix(0, n + n2, n + n2)
  V_red <- matrix(0, n + nnp12, n + nnp12)

  # Compute AUX
  AUX <- vec_1(nu_red + Psi_red %*% E_red)

  # Compute V and V_red if requested
  if (indic_compute_V == 1) {
    V_red <- AUX
    for (i in 1:100) {
      V_red <- AUX + Phi_tilde_red %*% V_red %*% t(Phi_tilde_red)
    }
    V <- M_ZY %*% V_red %*% t(M_ZY)
  }

  # Return results
  return(list(
    mu_tilde = mu_tilde,
    Phi_tilde = Phi_tilde,
    mu_tilde_red = mu_tilde_red,
    Phi_tilde_red = Phi_tilde_red,
    nu = nu,
    Psi = Psi,
    nu_red = nu_red,
    Psi_red = Psi_red,
    E = E,
    E_red = E_red,
    V = V,
    V_red = V_red,
    M_YZ = M_YZ,
    M_ZY = M_ZY
  ))
}



psi.GaussianQVAR <- function(u,psi.parameterization){
  # The model is:
  # x_t = mu + Phi·x_{t-1} + Sigma^(1/2)·eps_t,
  #     with eps_t ~ N(0,Id).
  # Hence, Sigma = Var_t(x_{t+1})
  # The extended state vector is: w_t = [x_t',vech(x_t·x_t')]'; if the
  #   dimension of x_t is n, the dimension of w_t is (n + n*(n+1)/2)

  mu <- matrix(psi.parameterization$mu,ncol=1)
  Phi <- psi.parameterization$Phi
  Sigma <- psi.parameterization$Sigma

  k <- dim(u)[2] # number of u's on which to evaluate the function
  n <- dim(Phi)[1]

  Id <- diag(n)

  a <- matrix(NaN,n+n*(n+1)/2,k)
  b <- matrix(NaN,k,1)

  Mx <- make_Mnx(n) # to go from vech to vec
  Kx <- make_Knx(n) # to go from vec to vech

  for(i in 1:k){
    v <- matrix(u[1:n,i],ncol=1)
    Vred <- u[(n+1):(n+n*(n+1)/2),i]
    vecV <- Mx %*% Vred
    V <- matrix(vecV,n,n)

    #Sigma_1_2V_1 <- solve(solve(Sigma) - 2*V)
    Sigma_1_2V_1 <- solve(Id - 2 * Sigma %*% V) %*% Sigma

    a1_i <- 2 * t(Phi) %*% V %*% Sigma_1_2V_1 %*% (v + 2*V %*% mu) +
      t(Phi) %*% v + 2 * t(Phi) %*% V %*% mu
    a2_i <- 2 * t(Phi) %*% V %*% Sigma_1_2V_1 %*% V %*% Phi +
      t(Phi) %*% V %*% Phi
    b_i <- .5 * t(v) %*% Sigma_1_2V_1 %*% v +
      2*t(v + V %*% mu) %*% Sigma_1_2V_1 %*% V %*% mu -
      .5*log(det(Id - 2*Sigma %*% V)) +
      t(v + V %*% mu) %*% mu

    a[,i] <- c(a1_i,Kx %*% c(a2_i))
    b[i]  <- b_i
  }

  return(list(a = a, b = b))
}


#det_complex <- function(M) prod(eigen(M, only.values = TRUE)$values)


make_Knx <- function(n) {
  S <- matrix(0, nrow = n, ncol = n)
  aux <- 1:(n * (n + 1) / 2)
  count <- 0

  # First loop to fill lower triangular part
  for (j in 1:n) {
    nb <- n - j + 1
    S[j:n, j] <- aux[(count + 1):(count + nb)]
    count <- count + nb
  }

  # Second loop to fill upper triangular part
  count <- 0
  for (j in 1:(n-1)) {
    nb <- n - j + 1
    S[j, (j + 1):n] <- aux[(count + 2):(count + nb)]
    count <- count + nb
  }

  # Create K matrix (note: original MATLAB code had M but should be K)
  K <- matrix(0, nrow = n * (n + 1) / 2, ncol = n^2)

  # Fill K matrix
  for (i in 1:(n^2)) {
    K[S[i], i] <- 1
  }

  return(K)
}

make_Mnx <- function(n) {
  S <- matrix(0, nrow = n, ncol = n)
  aux <- 1:(n * (n + 1) / 2)
  count <- 0

  # First loop to fill lower triangular part
  for (j in 1:n) {
    nb <- n - j + 1
    S[j:n, j] <- aux[(count + 1):(count + nb)]
    count <- count + nb
  }

  # Second loop to fill upper triangular part
  count <- 0
  for (j in 1:(n-1)) {
    nb <- n - j + 1
    S[j, (j + 1):n] <- aux[(count + 2):(count + nb)]
    count <- count + nb
  }

  # Create M matrix
  M <- matrix(0, nrow = n^2, ncol = n * (n + 1) / 2)

  # Fill M with 0.5 values
  for (i in 1:(n^2)) {
    M[i, S[i]] <- 0.5
  }

  # Handle diagonal elements
  diagS <- diag(S)
  G <- matrix(1:(n^2), nrow = n, ncol = n)
  diagG <- diag(G)

  for (i in 1:n) {
    M[diagG[i], diagS[i]] <- 1
  }

  return(M)
}

#
# mu <- matrix(c(0,1,.2),3,1)
# Phi <- .9*diag(3)
# Phi[2,1] <- .3
# aux <- matrix(.1*rnorm(9),3,3)
# Sigma <- aux %*% t(aux)
# #
# # res1 <- vcov_matrix_Quadratic(mu, Phi, Sigma, BigMats = make_BigMats(dim(Phi)[1]))
# # res2 <- make_matrices_cond_mean_variance_Quadratic(mu, Phi, Sigma, indic_compute_V=TRUE)
#
# model <- list(mu=mu,Phi=Phi,Sigma=Sigma)
# n <- dim(Phi)[1]
#
# Vred <- .1*rnorm(n*(n+1)/2)
#
# u <- matrix(c(.1*rnorm(n),Vred),ncol=1)
#
# psi_Quad <- psi.GaussianQVAR(u,psi.parameterization=model)
#
# model$n_w <- 9
# res <- compute_expect_variance(psi.GaussianQVAR,model,du = 10^(-6))
#
# res2 <- make_matrices_cond_mean_variance_Quadratic(model,indic_compute_V = TRUE)
#
# Model <- model
# Model$Sigma <- t(chol(model$Sigma))
# x <- simul.VAR(Model,nb.sim = 10000)
# xx <- (x %x% matrix(1,1,3)) * (matrix(1,1,3) %x% x)
# w <- cbind(x,xx)
# Sigma12 <- matrix(rnorm(4),2,2)
# aux <- matrix(rnorm(4),2,2)
# V <- aux + t(aux)
# Sigma12 %*% solve(diag(2) - 2*t(Sigma12) %*% V %*% Sigma12) %*% t(Sigma12)
# solve(solve(Sigma12%*%t(Sigma12)) - 2*V)
#
#
#
