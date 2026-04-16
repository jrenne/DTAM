# model <- (list(
#   alpha = matrix(c(.0,.0),ncol=1),
#   nu = matrix(c(0.2,0.2),ncol=1),
#   mu = matrix(c(.005,.005),ncol=1),
#   beta = NaN,
#   n_w = 2
# ))
# model$beta <- .95*diag(1/c(model$mu))
#
# set.seed(123)
#
# TT = 52*15
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
# W_save <- W
#
# plot(Nplus,type="l")
# lines(Nminus,col="red")
#
# xi0 <- .005
# xi1 <- matrix(c(0,0,.25/100,-.25/100),ncol=1)
#
# r <- xi0 + W %*% xi1
# plot(r,type="l",lwd=2,las=1)
# abline(h=0,col="grey",lty=3)
#
#
#
# #stop()
#
#
# # u <- matrix(rnorm(4*3),4,3)
# # psi.VARG_Poisson(u,psi.parameterization=model)
#
# model$n_w <- 4 # number of factors.
#
# HH <- 52*3
#
# res <- compute_AB_classical(xi0, xi1, kappa0 = NaN, kappa1 = NaN,
#                             H=HH, psi=psi.VARG_Poisson,
#                             psi.parameterization=model)
#
# yields <- matrix(1,TT,1) %*% matrix(res$b,nrow=1) +
#   W %*% matrix(res$a,dim(res$a)[1],dim(res$a)[3])
#
# t <- 530
# lines((t+1):(t+HH),yields[t,],col="red",lwd=2)
#
# F_affine <- compute_F_Shadow_affine(W=W,psi.VARG_Poisson,
#                                     psi.parameterization=model,
#                                     ell_bar=i_bar,b=xi0,a=xi1,c=0*xi1,
#                                     H=HH,
#                                     eps = 10^(-6), # to compute dG
#                                     max_x = 1000,
#                                     dx_statio = 10,
#                                     min_dx = 1e-06,
#                                     nb_x1 = 500)
#
# lines((t+1):(t+HH),cumsum(F_affine$F[t,])/(1:HH),col="blue",lwd=2)
#
