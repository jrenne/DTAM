


solve_learning <- function(model,max.iter=200){

  mu      <- as.matrix(model$mu,ncol=1)
  Phi     <- as.matrix(model$Phi)
  Sigma12 <- as.matrix(model$Sigma12)
  Omega12 <- as.matrix(model$Omega12)
  A       <- as.matrix(model$A)
  B       <- as.matrix(model$B,ncol=1)
  Lambda0 <- as.matrix(model$Lambda0)
  Lambda1 <- as.matrix(model$Lambda1)

  Sigma <- Sigma12 %*% t(Sigma12)
  Omega <- Omega12 %*% t(Omega12)

  m <- dim(A)[1]
  n <- dim(A)[2]

  Id_n <- diag(n)

  P0 <- Id_n
  P  <- P0
  old_P <- P

  for (i in 1:max.iter){
    Pstar <- Sigma + Phi %*% P %*% t(Phi)
    S   <- A %*% Pstar %*% t(A) + Omega
    S_1 <- solve(S)
    K   <- Pstar %*% t(A) %*% S_1
    P  <- (Id_n - Pstar %*% t(A) %*%
              solve(A %*% Pstar %*% t(A)) %*% A) %*% Pstar
    P.change <- P - old_P
    old_P <- P
  }

  R <- solve(Id_n + (Id_n - K %*% A) %*% Lambda0)

  model_sol    <- model
  model_sol$R  <- R
  model_sol$P <- P
  model_sol$K  <- K
  model_sol$S  <- S

  # Dynamics of (w_t',w_{t|t}')':

  AUX <- rbind(
    cbind(Id_n,Lambda0),
    cbind(0 * Id_n,Id_n))
  Gamma <- solve(AUX)

  mu_ww <- Gamma %*% rbind(mu,
                           R %*% mu)
  Phi_ww <- Gamma %*%
    rbind(cbind(Phi,Lambda1),
          cbind(R %*% K %*% A,
                R %*% (Id_n - K %*% A) %*% (Phi + Lambda1)))

  Sigma12_ww <- Gamma %*%
    rbind(cbind(Sigma12,matrix(0,n,m)),
          cbind(R %*% K %*% A %*% Sigma12,R %*% K %*% Omega12))

  Sigma_ww <- Sigma12_ww %*% t(Sigma12_ww)

  model_sol$mu_ww    <- mu_ww
  model_sol$Phi_ww   <- Phi_ww
  model_sol$Sigma_ww <- Sigma_ww
  model_sol$Sigma12_ww <- Sigma12_ww

  return(model_sol)
}


simul_model <- function(model_sol,H){

  mu_ww      <- model_sol$mu_ww
  Phi_ww     <- model_sol$Phi_ww
  Sigma12_ww <- t(chol(model_sol$Sigma_ww))

  n_ww <- length(mu_ww)
  n   <- n_ww/2

  x <- solve(diag(n_ww) - Phi_ww) %*% mu_ww
  X <- matrix(x,ncol=1)

  for(t in 2:H){
    x <- mu_ww + Phi_ww %*% x + Sigma12_ww %*% rnorm(dim(Sigma12_ww)[2])
    X <- cbind(X,x)
  }
  w_t  <- X[1:n,]
  w_tt <- X[(n+1):(2*n),]

  return(list(w_t  = w_t,
              w_tt = w_tt))
}

