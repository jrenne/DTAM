#' Top-down credit-model transforms and simulation
#'
#' The `psi.*TopDown` functions evaluate the one-step Laplace transform of the
#' latent loss/intensity state used in the top-down credit model.
#' `simul.TopDown()` simulates the state vector under the physical measure.
#'
#' @param u,u.y Matrix of transform arguments.
#' @param model List of top-down model parameters.
#' @param psi.y Optional function used to evaluate the transform of the `y`
#'   block.
#' @param nb_periods Number of simulated dates.
#' @param W0 Optional initial state.
#'
#' @return The `psi.*` functions return affine-transform coefficients. The
#'   simulator returns a list with unconditional moments (`EV`) and simulated
#'   paths for `y`, `n`, and the stacked state `w`.
#'
#' @examples
#' model <- list(
#'   alpha.n = c(0.2),
#'   beta.ny = matrix(0.1, 1, 1),
#'   beta.nn = matrix(0.3, 1, 1),
#'   beta.yy = matrix(0.2, 1, 1),
#'   beta.yn = matrix(0.05, 1, 1),
#'   alpha.y = c(0.1),
#'   nu.y = c(1),
#'   mu.y = c(0.4),
#'   alpha.m = matrix(0, 2, 1),
#'   n_w = 2
#' )
#' psi.y.TopDown(matrix(0.05, 1, 1), model)$b
#' sim <- simul.TopDown(model, nb_periods = 10)
#' dim(sim$all_w)
#'
#' @name TopDown_tools
#' @export
psi.w.TopDown <- function(u,model,psi.y=psi.y.TopDown){
  # This function evaluates the LT of the state vector w_t in the context of the
  # top-down approach of Gourieroux, Monfort, Mouabbi, and Renne.
  # "model" is a list containing the model parameters.
  # "psi.y" is the LT of y_{t+1} conditionally on w_t.

  alpha.n <- matrix(model$alpha.n,ncol=1)
  beta.ny <- model$beta.ny
  beta.nn <- model$beta.nn

  J   <- dim(alpha.n)[1] # number of segments in the economy.
  n.w <- dim(u)[1]
  n.y <- n.w - J

  u.y <- u[1:n.y,]
  u.n <- u[(n.y+1):(n.y+J),]

  ab_y <- psi.y.TopDown(u.y + t(beta.ny) %*% (exp(u.n) - 1),model)

  a.wy <- ab_y$a.yy
  a.wn <- ab_y$a.yn + t(beta.nn) %*% (exp(u.n) - 1)
  b.w  <- ab_y$b    + t(exp(u.n) - 1) %*% alpha.n

  return(list(b    = b.w,
              a.wy = a.wy,
              a.wn = a.wn,
              a    = rbind(a.wy,a.wn)))
}
#' @rdname TopDown_tools
#' @export
psiQ.w.TopDown <- function(u,model,psi.y=psi.y.TopDown){
  # Same as psiQ.w.TopDown, but under the risk-neutral (Q) measure

  k <- dim(u)[2]
  u_adjusted <- cbind(u + model$alpha.m %*% matrix(1,1,k),model$alpha.m)

  res <- psi.w.TopDown(u_adjusted,model,psi.y)
  b.w  <- res$b
  a.wy <- res$a.wy
  a.wn <- res$a.wn

  a.wy <- a.wy[,1:k] - a.wy[,k+1] %*% matrix(1,1,k)
  a.wn <- a.wn[,1:k] - a.wn[,k+1] %*% matrix(1,1,k)
  b.w  <- b.w[1:k,]  - b.w[k+1]

  return(list(b    = b.w,
              a.wy = a.wy,
              a.wn = a.wn,
              a    = rbind(a.wy,a.wn)))
}
#' @rdname TopDown_tools
#' @export
psi.y.TopDown <- function(u.y,model){
  # This function returns the LT of y_{t+1}|w_t.
  # It can be evaluated for different u.y's, i.e., u.y can be a matrix.
  # In this example, the conditional distribution of y_t is a compound
  #   Poisson-gamma (as in ARG processes).

  beta.yy <- model$beta.yy
  beta.yn <- model$beta.yn

  alpha.y <- model$alpha.y
  nu.y    <- model$nu.y
  mu.y    <- model$mu.y
  n.y <- dim(u.y)[1]
  k   <- dim(u.y)[2]
  mu_matrix <- matrix(mu.y,n.y,k)
  a.yy <- t(beta.yy) %*% ((u.y * mu_matrix)/(1 - u.y * mu_matrix))
  a.yn <- t(beta.yn) %*% ((u.y * mu_matrix)/(1 - u.y * mu_matrix))
  b.y <- t((u.y * mu_matrix)/(1 - u.y * mu_matrix)) %*% alpha.y -
    t(log(1 - u.y * mu_matrix)) %*% nu.y
  return(list(a.yy = a.yy,
              a.yn = a.yn,
              b    = b.y))
}

#' @rdname TopDown_tools
#' @export
simul.TopDown <- function(model,nb_periods,W0=NaN){
  # This function simulates the Top Down model. The conditional distribuion of
  # y_t is compound Poisson-Gamma (as in the ARG).
  beta.yy <- as.matrix(model$beta.yy)
  beta.yn <- as.matrix(model$beta.yn)
  beta.nn <- as.matrix(model$beta.nn)
  beta.ny <- as.matrix(model$beta.ny)

  alpha.n <- matrix(model$alpha.n,ncol=1)

  alpha.y <- matrix(model$alpha.y,ncol=1)
  nu.y    <- matrix(model$nu.y,ncol=1)
  mu.y    <- matrix(model$mu.y,ncol=1)

  n.y <- dim(alpha.y)[1]
  J   <- dim(beta.nn)[1]
  n.w <- n.y + J

  EV <- NULL
  if(is.na(W0)[1]){
    EV <- compute_expect_variance(psi.w.TopDown,model)
    Ew <- EV$Ew
    W0 <- matrix(Ew,ncol=1)
  }
  W <- matrix(W0,ncol=1)
  y <- matrix(W0[1:n.y],ncol=1)
  n <- matrix(W0[(n.y+1):n.w],ncol=1)
  all_y <- y
  all_n <- n

  for(t in 2:nb_periods){
    # y part:
    param.pois <- rpois(n.y, alpha.y + beta.yy %*% y + beta.yn %*% n)
    y <- rgamma(n.y,shape = param.pois + nu.y, scale = mu.y)
    all_y <- cbind(all_y,y)
    # n part:
    n <- rpois(J, alpha.n + beta.ny %*% y + beta.nn %*% n)
    all_n <- cbind(all_n,n)
  }

  return(list(
    EV = EV,
    all_y = all_y,
    all_n = all_n,
    all_w = rbind(all_y,
                  all_n)
  ))
}


compute_H_bar_TopDown <- function(model,gamma,H=10,W=NaN,
                                  indic_Q = TRUE){
  # gamma has to be a vector with the same dimension as w_t
  # model contains the parameterization of the VARG model.
  # Prices are computed for all segments.

  gamma <- matrix(gamma,ncol=1)
  xi1 <- matrix(model$xi1,ncol=1)
  xi0 <- model$xi0

  epsilon <- 10^(-7) # used to numerically compute gradient

  # Build u1 and u2 -------------
  # (1st column: for H_upper; 2nf column: for H_lower,
  #  3rd and 4th columns: are the same, for H0)
  u1 <- cbind(gamma*epsilon,
              gamma*0,
              gamma*0,
              gamma*0)
  u2 <- cbind(gamma*epsilon - xi1,
              gamma*epsilon - xi1,
              gamma*0       - xi1,
              gamma*0       - xi1)

  if(indic_Q){
    psi <- psiQ.w.TopDown
  }else{
    psi <- psi.w.TopDown
  }
  varphi <- reverse.MHLT(psi,
                         u1 = u1,
                         u2 = u2,
                         H = H,model)
  nw <- dim(varphi$A)[1]

  da_du <- (varphi$A[,1:2,] - varphi$A[,3:4,])/epsilon
  db_du <- (varphi$B[,1:2,] - varphi$B[,3:4,])/epsilon
  db_du <- array(db_du,c(1,2,H))

  a_exp <- varphi$A[,1:2,] - array(xi1,c(nw,2,H))
  b_exp <- array(varphi$B[,1:2,],c(1,2,H)) -
    array(xi0*(1:H)%x%c(1,1),c(1,2,H))

  # Risk-free bond prices:
  a_exp_rf <- varphi$A[,3:4,] - array(xi1,c(nw,2,H))
  b_exp_rf <- array(varphi$B[,3:4,],c(1,2,H)) -
    array(xi0*(1:H)%x%c(1,1),c(1,2,H))
  a_exp_rf <- a_exp_rf[,1,]
  b_exp_rf <- b_exp_rf[,1,]

  if(!is.na(W[1])){
    if(dim(W)[2]!=nw){
      # to make sure W is T x nw
      W <- t(W)
    }
    T <- dim(W)[1]
    H_upper <- matrix(0,T,H)
    H_lower <- matrix(0,T,H)
    B       <- matrix(0,T,H)
    for(h in 1:H){
      H_upper[,h] <- exp(b_exp[1,1,h] + W %*% a_exp[,1,h]) *
        (db_du[1,1,h] + W %*% da_du[,1,h])
      H_lower[,h] <- exp(b_exp[1,2,h] + W %*% a_exp[,2,h]) *
        (db_du[1,2,h] + W %*% da_du[,2,h])
      B[,h] <- exp(b_exp_rf[h] + W %*% a_exp_rf[,h])
    }
  }else{
    H_upper <- NaN
    H_lower <- NaN
    B       <- NaN
  }
  return(list(da_du = da_du,
              db_du = db_du,
              a_exp = a_exp,
              b_exp = b_exp,
              a_exp_rf = a_exp_rf,
              b_exp_rf = b_exp_rf,
              H_upper = H_upper,
              H_lower = H_lower,
              B = B))
}


#' Top-down credit pricing objects
#'
#' `compute_CDS_TopDown()` prices index-segment CDS contracts in the top-down
#' model. `compute_G()` evaluates the truncated-payoff building block used for
#' loss-function transforms. `compute_tranche()` combines these ingredients into
#' tranche prices.
#'
#' @param model List of top-down model parameters.
#' @param H Maximum horizon.
#' @param W Matrix of state-vector values.
#' @param RR Recovery rate.
#' @param indic_Q Logical. If `TRUE`, work under the risk-neutral transform.
#' @param gamma,v Inputs defining the truncated payoff.
#' @param b.matrix Threshold matrix, one row per maturity.
#' @param indic_upper Logical determining whether the upper or lower transform is
#'   computed.
#' @param max_x,dx_statio,min_dx,nb_x1 Numerical integration settings.
#' @param j Segment index for tranche pricing.
#' @param vector.of.a Attachment/detachment points.
#'
#' @return `compute_CDS_TopDown()` returns a `T x J x H` array of CDS spreads.
#'   `compute_G()` returns the array of truncated-payoff values. 
#'   `compute_tranche()` returns a list with `tranche_prices`, `numerator`, and
#'   `denominat`.
#'
#' @examples
#' model <- list(
#'   alpha.n = c(0.2),
#'   beta.ny = matrix(0.1, 1, 1),
#'   beta.nn = matrix(0.3, 1, 1),
#'   beta.yy = matrix(0.2, 1, 1),
#'   beta.yn = matrix(0.05, 1, 1),
#'   alpha.y = c(0.1),
#'   nu.y = c(1),
#'   mu.y = c(0.4),
#'   alpha.m = matrix(0, 2, 1),
#'   xi0 = 0,
#'   xi1 = matrix(0, 2, 1),
#'   I = 1,
#'   n_w = 2
#' )
#' W <- matrix(c(0.3, 0), 1, 2)
#' cds <- compute_CDS_TopDown(model, H = 2, W = W)
#' dim(cds)
#'
#' @name TopDown_pricing
compute_CDS_TopDown <- function(model,
                                H=10,
                                W,
                                RR=.5,
                                indic_Q=TRUE){

  alpha.n <- matrix(model$alpha.n,ncol=1)

  J   <- length(model$alpha.n) # number of segments in the economy.
  n.y <- length(model$alpha.y) # number of y variables.
  n.w <- J + n.y

  if(dim(W)[2]!=n.w){# to make sure W is of dimension T x n.w
    W <- t(W)
  }
  T <- dim(W)[1]

  CDS <- array(0,c(T,J,H))

  for(j in 1:J){
    I_j <- model$I[j]
    ej_tilde <- matrix(0,n.w,1)
    ej_tilde[n.y+j] <- 1
    res <- compute_H_bar_TopDown(model,ej_tilde,H=H,W,
                                 indic_Q = indic_Q)
    numerator <- 0
    denominat <- 0
    for(h in 1:H){
      numerator <- numerator + res$H_upper[,h] - res$H_lower[,h]
      denominat <- denominat + I_j * res$B[,h] - res$H_lower[,h]
      CDS[,j,h] <- (1-RR)*(numerator/denominat)
    }
  }

  return(CDS)
}


varphi4G_TopDown <- function(x,parameterization,H){

  model   <- parameterization$model
  gamma   <- parameterization$gamma
  indic_Q     <- parameterization$indic_Q # if FALSE, use physical LT, otherwise Q
  indic_upper <- parameterization$indic_upper # if TRUE, compute G_upper, otherwise G_lower

  v <- parameterization$v # determine thresholds (dimension nw x 1)
  i.v.x <- matrix(v,ncol=1) %*% matrix(c(1i*x),nrow=1)

  nw <- dim(i.v.x)[1]
  q  <- dim(i.v.x)[2]

  u2 <- matrix(gamma - model$xi1, nw, q) + i.v.x
  u1 <- i.v.x
  if(indic_upper){
    u1 <- u1 + matrix(gamma, nw, q)
  }

  if(indic_Q){
    psi <- psiQ.w.TopDown
  }else{
    psi <- psi.w.TopDown
  }

  res_reverse <- reverse.MHLT(psi,
                              u1 = u1,
                              u2 = u2,
                              H = H,
                              psi.parameterization = model)

  A <- res_reverse$A
  A <- A - array(matrix(1,q*H,1) %x% matrix(model$xi1,ncol=1),c(nw,q,H))
  B <- res_reverse$B
  B <- B - model$xi0*array(matrix((1:H),ncol=1) %x% matrix(1,q,1),c(1,q,H))

  return(list(A = res_reverse$A,
              B = matrix(res_reverse$B,nrow=1)))
}


truncated.payoff <- function(W, # values of w_t
                             b.matrix, # thresholds (H x k, 1 row per maturity)
                             varphi,
                             parameterization,
                             max_x = 500,
                             dx_statio = .1,
                             min_dx = 1e-05,
                             nb_x1 = 1000){
  # This function computes:
  #       E[exp(u0·w_t+...+uh·w_{t+h})·1_{v0·w_t+...+vh·w_{t+h} < b}],
  # for different horizons and values of b.
  # b.matrix is organized as follows: H * k; the h^th row gives the b values for
  #       horizon h.
  # "W" is of dimension T x n.w; it contains values of the state vector w
  #       at which we want to evaluate the prices.
  # "varphi" is a function that takes three arguments:
  #       (i) u and
  #       (ii) parameterization (that contains at least "model"),
  #       (iii) H (max maturity);
  # it returns matrices A(u) and B(u) that are such that the "price" of
  #       the payoff is given by exp(A(u)'w + B(u)),
  #       where A(u) is of dimension n.u x n.w,
  #       and B(u) is of dimension n.u x 1
  # NOTE: It is varphi that determines the ui's and vi's.

  if(is.null(parameterization$model$n_w)){
    print("Error: 'parameterization' should include 'model', itself including 'n_w' (size of state vector).")
    return(1)
  }

  n_w <- parameterization$model$n_w

  dx1 <- exp(seq(log(min_dx),log(dx_statio),length.out=nb_x1))
  max_x1 <- tail(cumsum(dx1),1)
  nb_x2 <- (max_x - max_x1)/dx_statio
  dx <- c(dx1,rep(dx_statio,nb_x2))
  x <- matrix(cumsum(dx),ncol=1)

  nb_x <- length(x)

  H <- dim(b.matrix)[1]

  res_varphi  <- varphi(x,parameterization,H)
  res_varphi0 <- varphi(0,parameterization,H)

  # Adjust for current short-term interest rate:
  B <- matrix(res_varphi$B,1,nb_x*H)
  A <- matrix(res_varphi$A,n_w,nb_x*H)
  B0 <- matrix(res_varphi0$B,nrow=1)
  A0 <- matrix(res_varphi0$A,n_w,H)

  if(dim(W)[2]!=n_w){# to make sure W is of dimension T x n.w
    W <- t(W)
  }
  T <- dim(W)[1]

  psi_eval  <- exp(t(A)  %*% t(W) + matrix(B, nb_x*H,T))
  psi_eval0 <- Re(exp(t(A0) %*% t(W) + matrix(B0,H,T)))

  # modify format to have them nb_x x (T*H):
  Psi_eval  <- matrix(NaN,nb_x,T*H)
  for(h in 1:H){
    Psi_eval[,(T*(h-1)+1):(T*h)]  <- psi_eval[(nb_x*(h-1)+1):(nb_x*h),]
  }

  dx <- matrix(x - c(0,x[1:(nb_x-1)]),nb_x,T*H)
  x_mid <- matrix(x,nb_x,T*H)

  all_prices <- matrix(0,T*H,dim(b.matrix)[2])

  count_b <- 0
  for(b in 1:dim(b.matrix)[2]){
    count_b <- count_b + 1
    aux <- Im(Psi_eval *
                exp(-1i*matrix(x,ncol=1) %*%
                      (matrix(b.matrix[,count_b],nrow=1) %x% matrix(1,1,T))))
    fx <- aux/x_mid * dx
    f  <- 1/2*c(t(psi_eval0)) - 1/pi * apply(fx,2,sum)
    all_prices[,count_b] <- f
  }

  All_prices <- array(NaN,c(T,dim(b.matrix)[2],H))
  for(h in 1:H){
    All_prices[,,h] <- all_prices[(T*(h-1)+1):(T*h),]
  }

  return(All_prices)
}


#' @rdname TopDown_pricing
#' @export
compute_G <- function(model,W,
                      gamma,v,
                      b.matrix, # attachment/detachment points
                      H, # maturity
                      indic_Q = TRUE, # computed under Q.
                      indic_upper = TRUE, # otherwise, computed G_lower
                      max_x = 2000,
                      dx_statio = 2,
                      min_dx = 10^(-6),
                      nb_x1=1000){
  # This computes functions G_lower and G_upper

  parameterization <- list(model=model,
                           indic_Q = indic_Q,
                           indic_upper = indic_upper,
                           gamma = gamma,
                           v = v)

  # P <- truncated.payoff(W,
  #                       v,vector.of.b,
  #                       H,
  #                       varphi = varphi4G_TopDown,
  #                       parameterization,
  #                       max_x = max_x,
  #                       dx_statio = dx_statio,
  #                       min_dx = min_dx,
  #                       nb_x1=nb_x1)
  P <- truncated.payoff(W,
                        b.matrix,
                        varphi = varphi4G_TopDown,
                        parameterization,
                        max_x = max_x,
                        dx_statio = dx_statio,
                        min_dx = min_dx,
                        nb_x1=nb_x1)

  return(P)
}

#' @rdname TopDown_pricing
#' @export
compute_tranche <- function(model,W,
                            j,H,
                            vector.of.a, # attachment/detachment points
                            RR = 0.5, # recovery rate
                            max_x = 500,
                            dx_statio = .1,
                            min_dx = 1e-05,
                            nb_x1 = 1000,
                            indic_Q = TRUE # if FALSE; then everyting computed under P
){
  # This function computes the price of credit tranche products in
  #     the context of the top-down credit approach.
  # It does that for one segment only (segment j).
  # vector.of.a contains the attachment/detachment points. For instance,
  #     if vector.of.a =[0,.03,.06], then we consider 2 tranches:
  #     [0,3%] and [3%,6%].
  # W is a matrix of dimension T x n.w containing the consider values
  #     of the state vector.

  n_w <- model$n_w
  if(dim(W)[2]!=n_w){# to make sure W is of dimension T x n.w
    W <- t(W)
  }
  T <- dim(W)[1]

  vector.of.a.bar <- vector.of.a * model$I[j] / (1 - RR)
  nb_a <- length(vector.of.a)

  # Prepare arrays of a.bar and b.bar
  # b.bar:
  aux <- vector.of.a.bar[2:nb_a]
  aux <- matrix(1,T,1) %*% matrix(aux,1,nb_a-1)
  array.b.bar <- array(aux,c(T,nb_a-1,H))
  # a.bar:
  aux <- vector.of.a.bar[1:(nb_a-1)]
  aux <- matrix(1,T,1) %*% matrix(aux,1,nb_a-1)
  array.a.bar <- array(aux,c(T,nb_a-1,H))

  n.w <- model$n_w
  J <- dim(as.matrix(model$beta.nn))[1]
  n.y <- n.w - J

  ej_tilde <- matrix(0,n.w,1)
  ej_tilde[n.y + j] <- 1

  G_upper_0_ej <- compute_G(model,
                            W=W,
                            gamma = 0*ej_tilde,
                            v = ej_tilde,
                            b.matrix=t(matrix(vector.of.a.bar,
                                              length(vector.of.a.bar),H)), # attachment/detachment points
                            H=H,
                            indic_Q = indic_Q,
                            indic_upper = TRUE,
                            max_x = max_x,dx_statio = dx_statio,
                            min_dx = min_dx,nb_x1 = nb_x1)
  G_lower_0_ej <- compute_G(model,
                            W=W,
                            gamma = 0*ej_tilde,
                            v = ej_tilde,
                            b.matrix=t(matrix(vector.of.a.bar,
                                              length(vector.of.a.bar),H)), # attachment/detachment points
                            H=H,
                            indic_Q = indic_Q,
                            indic_upper = FALSE,
                            max_x = max_x,dx_statio = dx_statio,
                            min_dx = min_dx,nb_x1 = nb_x1)

  u <- 10^(-4)
  G_upper_uej_ej <- compute_G(model,
                              W=W,
                              gamma = u*ej_tilde,
                              v = ej_tilde,
                              b.matrix=t(matrix(vector.of.a.bar,
                                                length(vector.of.a.bar),H)), # attachment/detachment points
                              H=H,
                              indic_Q = indic_Q,
                              indic_upper = TRUE,
                              max_x = max_x,dx_statio = dx_statio,
                              min_dx = min_dx,nb_x1 = nb_x1)
  G_lower_uej_ej <- compute_G(model,
                              W=W,
                              gamma = u*ej_tilde,
                              v = ej_tilde,
                              b.matrix=t(matrix(vector.of.a.bar,
                                                length(vector.of.a.bar),H)), # attachment/detachment points
                              H=H,
                              indic_Q = indic_Q,
                              indic_upper = FALSE,
                              max_x = max_x,dx_statio = dx_statio,
                              min_dx = min_dx,nb_x1 = nb_x1)

  dG_upper <- (G_upper_uej_ej - G_upper_0_ej)/u
  dG_lower <- (G_lower_uej_ej - G_lower_0_ej)/u

  G_upper_0_ej_a <- G_upper_0_ej[,1:(nb_a-1),]
  G_upper_0_ej_b <- G_upper_0_ej[,2:nb_a,]

  dG_upper_a <- dG_upper[,1:(nb_a-1),]
  dG_upper_b <- dG_upper[,2:nb_a,]

  dG_lower_a <- dG_lower[,1:(nb_a-1),]
  dG_lower_b <- dG_lower[,2:nb_a,]

  numerator <- dG_upper_b - dG_upper_a - (dG_lower_b - dG_lower_a)
  denominat <- dG_upper_a - array.a.bar*G_upper_0_ej_a +
    array.b.bar*G_upper_0_ej_b - dG_upper_b

  for(h in 2:H){
    numerator[,,h] <- numerator[,,h-1] + numerator[,,h]
    denominat[,,h] <- denominat[,,h-1] + denominat[,,h]
  }

  tranche_prices <- numerator/denominat

  return(list(tranche_prices=tranche_prices,
              numerator=numerator,
              denominat=denominat))
}


compute_affine_payoff_TopDown <- function(model,gamma,H,
                                          indic_Q = TRUE,
                                          W=NaN){
  # This functions computes the date-t price of exp(gamma'w_{t+h}) in
  # the context of the top-down approach.

  n_w <- model$n_w

  u2 <- matrix(- model$xi1, ncol=1)
  u1 <- matrix(gamma,n_w,1)

  if(indic_Q){
    psi <- psiQ.w.TopDown
  }else{
    psi <- psi.w.TopDown
  }

  res_reverse <- reverse.MHLT(psi,
                              u1 = u1,
                              u2 = u2,
                              H = H,
                              psi.parameterization = model)

  A <- matrix(res_reverse$A,n_w,H)
  B <- matrix(res_reverse$B,1,H)

  A <- A - matrix(model$xi1,ncol=1) %*% matrix(1,1,H)
  B <- B - matrix(model$xi0,ncol=1) %*% matrix(1:H,1,H)

  a <- A/t(matrix(1:H,H,n_w))
  b <- B/t(matrix(1:H,H,1))

  if(!is.na(W[1])){
    # make sure W of the size T x n_w
    if(dim(W)[1]==n_w){
      W <- t(W)
    }
    T <- dim(W)[1]
    P <- exp(W %*% A + matrix(1,T,1) %*% B)
    r <- -log(P)/t(matrix(1:H,H,T))
  }else{
    P <- NaN
    r <- NaN
  }

  return(list(A = A,
              B = B,
              a = a,
              b = b,
              P = P,
              r = r))
}

#' Rating-migration bond pricing and simulation
#'
#' `compute_bond_price_Ratings()` prices risky and risk-free bonds in a model
#' with stochastic rating migration. `simul.rating.migration()` simulates a
#' rating path conditional on an exogenous state trajectory.
#'
#' @param model List of model parameters.
#' @param Y Matrix of factor realizations.
#' @param H Maximum maturity.
#' @param psi One-step Laplace transform of the factor process.
#' @param tau_ini Initial rating for the simulated migration path.
#'
#' @return `compute_bond_price_Ratings()` returns a list containing risk-free and
#'   risky bond prices, yields, and spreads. `simul.rating.migration()` returns
#'   a vector of simulated ratings.
#'
#' @examples
#' model <- list(
#'   kappa0 = matrix(c(0.01, 0), 1, 2),
#'   kappa1 = matrix(c(0.1, 0), 1, 2),
#'   xi0 = 0.01,
#'   xi1 = matrix(0.1, 1, 1),
#'   V = diag(2),
#'   psi.parameterization = list(mu = matrix(0, 1, 1),
#'                               Phi = matrix(0.8, 1, 1),
#'                               Sigma12 = matrix(0.1, 1, 1))
#' )
#' Y <- matrix(c(0, 0.1, -0.1), 3, 1)
#' res <- compute_bond_price_Ratings(model, Y, H = 2, psi = psi.GaussianVAR)
#' dim(res$RF_bond_price)
#' simul.rating.migration(model, Y, tau_ini = 1)
#'
#' @name Ratings_tools
#' @export
