#
# rho <- .98
# mu <- .16
#
# mu <- matrix( mu*(1-rho),2,1)
# Phi <- rho*diag(2)
# Sigma12 <- .02*matrix(c(1,-1),2,1)
# Sigma <- Sigma12 %*% t(Sigma12)
#
# model <- list(mu=mu,Phi=Phi,Sigma=Sigma,Sigma12=Sigma12,n_w=2)
#
# TT = 400
#
# simX <- simul.GVAR(model,nb.sim = TT)
#
# par(mfrow=c(2,2))
# par(plt=c(.1,.95,.1,.9))
# plot(simX[,1],type="l",ylim=c(min(0,min(simX)),max(simX)))
# lines(simX[,2],type="l",col="red")
# abline(h=0)
#
# lambda_plus  <- simX[,1]^2
# lambda_minus <- simX[,2]^2
#
# plot(lambda_plus,type="l",ylim=c(0,max(simX^2)))
# lines(lambda_minus,type="l",col="red")
# abline(h=0)
#
# DeltaNplus <- rpois(n = TT,lambda_plus)
# DeltaNminus <- rpois(n = TT,lambda_plus)
#
# Nplus  <- cumsum(DeltaNplus)
# Nminus <- cumsum(DeltaNminus)
#
# W <- cbind(simX,
#            ((simX %x% matrix(1,1,2))*(matrix(1,1,2) %x% simX))[,c(1,2,4)],
#            Nplus,Nminus)
#
#
# plot(Nplus,type="l")
# lines(Nminus,col="red")
#
# r <- 1 + .25*(Nplus - Nminus)
# plot(r,type="l",lwd=2)
# abline(h=0)
#
#
# n <- dim(model$Phi)[1]
#
# xi0 <- .01
# xi1 <- matrix(c(rep(0,n),
#                 rep(0,n*(n+1)/2),
#                 .25/100,-.25/100),ncol=1)
# H <- 20
#
# model$Omega_plus  <- diag(c(1,0))
# model$Omega_minus <- diag(c(0,1))
#
# u <- - xi1
# psi.QPoisson(u,model)
#
# model$n_w <- n + n*(n+1)/2 + 2
#
# res <- compute_AB_classical(xi0, xi1, kappa0 = NaN, kappa1 = NaN,
#                             H=H, psi=psi.QPoisson,
#                             psi.parameterization=model)
#
# yields <- matrix(1,TT,1) %*% matrix(res$b,nrow=1) +
#   W %*% matrix(res$a,dim(res$a)[1],dim(res$a)[3])
#
# t <- 40
# lines((t+1):(t+H),100*yields[t,],col="red",lwd=2)
#
#
# parameterization <- list(
#   model = model,
#   u = matrix(0*xi1,ncol=1),
#   v = matrix(1*xi1,ncol=1)
# )
#
# # varphi4G_SR_QPoisson(1,parameterization,H)
# # res_truncated <- truncated.payoff(W,b.matrix = matrix(-xi0,H,1),
# #                                   varphi = varphi4G_SR_QPoisson,
# #                                   parameterization = parameterization,
# #                                   max_x = 1000,
# #                                   dx_statio = 1,
# #                                   min_dx = 1e-07,
# #                                   nb_x1 = 2000)
# # h <- 10
# # plot(res_truncated[,1,h],type="l")
#
# #stop()
#
# i_bar <- 0
#
# rho <- .98
# mu <- .02
# mu <- matrix( mu*(1-rho),2,1)
# Phi <- rho*diag(2)
# Sigma12 <- .015*matrix(c(1,-1),2,1)
# Sigma <- Sigma12 %*% t(Sigma12)
#
# model <- list(mu=mu,Phi=Phi,Sigma=Sigma,Sigma12=Sigma12,n_w=2)
#
# TT = 400
# simX <- simul.GVAR(model,nb.sim = TT)
#
# par(mfrow=c(2,1))
# abline(h=0,col="grey")
# plot(simX[,1],type="l")
#
# xi1 <- matrix(c(1,0),ncol=1)
# xi0 <- 0
#
# parameterization <- list(
#   model = model,
#   u = matrix(0*xi1,ncol=1),
#   v = matrix(1*xi1,ncol=1)
# )
# res_truncated0 <- truncated.payoff(simX,b.matrix = matrix(i_bar-xi0,H,1),
#                                    varphi = varphi4G_SR_Gaussian,
#                                    parameterization = parameterization,
#                                    max_x = 2000,
#                                    dx_statio = 10,
#                                    min_dx = 1e-06,
#                                    nb_x1 = 1000)
# eps <- 10^(-6)
# parameterization$u <- matrix(eps*xi1,ncol=1)
# res_truncatedeps <- truncated.payoff(simX,b.matrix = matrix(i_bar-xi0,H,1),
#                                      varphi = varphi4G_SR_Gaussian,
#                                      parameterization = parameterization,
#                                      max_x = 2000,
#                                      dx_statio = 10,
#                                      min_dx = 1e-06,
#                                      nb_x1 = 1000)
# G0   <- matrix(res_truncated0,nrow=TT)
# Geps <- matrix(res_truncatedeps,nrow=TT)
#
# dG <- (Geps - G0)/eps
#
# # h <- 10
# # plot(G0[,h],type="l",ylim=c(0,1))
# # lines(G0[,1],col="red")
#
# model$n_w <- 2
# res_EV <- compute_expect_variance_H(compute_expect_variance(psi.GaussianVAR,model),
#                                     theta1 = t(xi1),theta2 = t(0*xi1),H = H)
#
# C0 <- matrix(res_EV$all_C0,nrow=1)
# C1 <- matrix(res_EV$all_C1,ncol=H)
# E_aW <- matrix(1,TT,1) %*% C0 + simX %*% C1
#
# res_EV2 <- compute_expect_variance_H(compute_expect_variance(psi.GaussianVAR,model),
#                                      theta1 = t(xi1),theta2 = t(xi1),H = H)
# delta_sigma2 <- c(res_EV2$all_Gamma0) - c(0,c(res_EV2$all_Gamma0)[1:(H-1)])
#
# F <- xi0*(1 - G0) + E_aW - dG - 1/2*(1-G0)*t(matrix(delta_sigma2,H,TT))
#
# h <- 20
#
# plot(F[,h],type="l")
#
# # Compare with Extended Wu-Xia:
# b <- xi0
# a <- xi1
# c <- 0*xi1
# res_GaussianLemma <- apply_Lemma(model$mu,
#                                  model$Phi,
#                                  model$Sigma12,i_bar=i_bar,b,a,c,simX,max_h=H)
# lines(res_GaussianLemma$F_c_equal_0[,h],col="red")
#
# #plot(simX[,1],type="l")
# plot(1-res_GaussianLemma$ptn[,h],type="l")
# lines(G0[,h],col="red")
#
# sigma <- t(sqrt(matrix(res_EV$all_Gamma0,H,TT)))
# mu <- E_aW
# beta <- - mu/sigma
# EaW1_aW <- pnorm(beta)*mu - dnorm(beta)*sigma
#
# plot(EaW1_aW[,h],type="l") # Exact formula
# lines(dG[,h],col="red") # Numerical approx
#
#
#
# g <- function(x){
#   return(x * pnorm(x) + dnorm(x))
# }
#
# all_sigma_n <- res_GaussianLemma$all_sigma_n
# all_b_n     <- res_GaussianLemma$all_b_n
# all_b_bar_n <- res_GaussianLemma$all_b_bar_n
# all_a_n     <- res_GaussianLemma$all_a_n
#
# F_lemma <- (matrix(1,TT,1) %*% t(all_sigma_n)) * g(
#   (matrix(1,TT,1) %*% t(all_b_n) +
#      simX %*% t(all_a_n))/
#     (matrix(1,TT,1) %*% t(all_sigma_n))
# )
#
# F_lemma_bis <- (matrix(1,TT,1) %*% t(all_sigma_n)) * g(
#   (matrix(1,TT,1) %*% t(all_b_bar_n) +
#      simX %*% t(all_a_n))/
#     (matrix(1,TT,1) %*% t(all_sigma_n))) +
#   (matrix(1,TT,1) %*% t(all_b_n - all_b_bar_n)) * pnorm(
#     (matrix(1,TT,1) %*% t(all_b_bar_n) +
#        simX %*% t(all_a_n))/
#       (matrix(1,TT,1) %*% t(all_sigma_n)))
#
#
# plot(F[,h],type="l")
# lines(F_lemma[,h],col="red")
# lines(F_lemma_bis[,h],col="blue")
#
# F_affine <- compute_F_Shadow_affine(W=simX,psi.GaussianVAR,
#                                     psi.parameterization=model,
#                                     ell_bar=i_bar,b=xi0,a=xi1,c=0*xi1,
#                                     H,
#                                     eps = 10^(-6), # to compute dG
#                                     max_x = 2000,
#                                     dx_statio = 10,
#                                     min_dx = 1e-06,
#                                     nb_x1 = 1000)
# lines(F_affine$F[,h],col="green")
#
# # plot(F_affine$E_cW[,10])
# # lines(E_aW[,10])
#
#
# stop()

# model <- (list(
#   alpha = matrix(c(.0,.0),ncol=1),
#   nu = matrix(c(.1,.1),ncol=1),
#   mu = matrix(c(.02,.02),ncol=1),
#   beta = NaN,
#   n_w = 2
# ))
# model$beta <- .95*diag(1/c(model$mu))
#
# set.seed(123)
#
# TT = 480
#
# simX <- t(simul_VARG(model,nb_periods = TT-1))
#
# par(mfrow=c(2,2))
#
# plot(simX[,1],type="l")
# lines(simX[,2],col="red")
#
# DeltaNplus <- rpois(n = TT,simX[,1])
# DeltaNminus <- rpois(n = TT,simX[,2])
#
# Nplus  <- cumsum(DeltaNplus)
# Nminus <- cumsum(DeltaNminus)
#
# W <- cbind(simX,
#            Nplus,Nminus)
#
# plot(Nplus,type="l")
# lines(Nminus,col="red")
#
# xi0 <- .01
# xi1 <- matrix(c(0,0,.25/100,-.25/100),ncol=1)
#
# r <- xi0 + W %*% xi1
# plot(r,type="l",lwd=2,las=1)
# abline(h=0,col="grey",lty=3)
#
# # u <- matrix(rnorm(4*3),4,3)
# # psi.VARG_Poisson(u,psi.parameterization=model)
#
# model$n_w <- 4 # number of factors.
#
# H <- 100
#
# res <- compute_AB_classical(xi0, xi1, kappa0 = NaN, kappa1 = NaN,
#                             H=H, psi=psi.VARG_Poisson,
#                             psi.parameterization=model)
#
# yields <- matrix(1,TT,1) %*% matrix(res$b,nrow=1) +
#   W %*% matrix(res$a,dim(res$a)[1],dim(res$a)[3])
#
#
# HH <- 100
# F_affine <- compute_F_Shadow_affine(W=W,psi.VARG_Poisson,
#                                     psi.parameterization=model,
#                                     ell_bar=0,b=xi0,a=xi1,c=0*xi1,
#                                     H=HH,
#                                     eps = 10^(-6), # to compute dG
#                                     max_x = 2000,
#                                     dx_statio = 10,
#                                     min_dx = 1e-06,
#                                     nb_x1 = 1000)
#
# t <- 20
# lines((t+1):(t+H),yields[t,],col="red",lwd=2)
# lines((t+1):(t+H),F_affine$F[t,],col="blue")
#
# #plot(F_affine$delta_sigma2_c_a[,100],type="l")
