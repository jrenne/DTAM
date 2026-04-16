#
#
# psi.stocks2fact <- function(u, psi.parameterization){
#   alpha <- psi.parameterization$alpha
#   beta  <- psi.parameterization$beta
#   phi   <- psi.parameterization$phi
#   sigma <- psi.parameterization$sigma
#   chi   <- psi.parameterization$chi
#
#   u1 <- c(u[1,])
#   u2 <- c(u[2,])
#
#   a <- rbind(u1*phi,
#              chi*u1^2/2 + u2*mu*beta/(1-mu*u2))
#
#   b <- matrix(1/2*u1^2*sigma^2 + u2*mu*alpha/(1-mu*u2),ncol=1)
#
#   return(list(
#     a = a, b = b
#   ))
# }
#
# simul.w <- function(model, nb.sim, w1_ini = 0, w2_ini = 0){
#
#   all_w2 <- simul.ARG(nb.sim = nb.sim,
#                       mu = model$mu,rho = model$mu*model$beta, nu = 0,
#                       alpha = model$alpha)
#
#   all_w1 <- NULL
#   w1 <- w1_ini
#   w2 <- w2_ini
#
#
#   for(t in 1:nb.sim){
#     w1 <- model$phi * w1 + sqrt(model$sigma^2 + model$chi*w2)*rnorm(1)
#     w2 <- all_w2[t]
#
#     all_w1 <- c(all_w1,w1)
#   }
#
#   return(cbind(all_w1,all_w2))
# }
#
#
# freq <- 4 # quarterly data
#
# alpha <- .05
# rho   <- .9
# mu    <- 2
# beta  <- rho/mu
# sigma <- .002
# chi   <- .00005/freq
# phi   <- .8
# kappa <- .000
# mu_pi0 <- .03/freq
#
# model <- list(
#   alpha = alpha,
#   mu    = mu,
#   beta  = beta,
#   sigma = sigma,
#   chi   = chi,
#   phi   = phi,
#   kappa = kappa,
#   mu_pi0 = mu_pi0,
#   n_w = 2
# )
#
# set.seed(2)
#
# nb.sim <- 200
# sim.w <- simul.w(model, nb.sim = nb.sim, w1_ini = 0, w2_ini = 0)
#
# w1 <- sim.w[,1]
# w2 <- sim.w[,2]
#
# pi_t <- model$mu_pi0 + w1 - model$kappa * w2
# annualized_pi_t <- freq * pi_t
#
# parameterization = list()
# parameterization$model <- model
# parameterization$xi0 = 0
# parameterization$xi1 = matrix(0,2,1)
# parameterization$mu_pi0 = model$mu_pi0
# parameterization$mu_pi1 = matrix(c(1,-model$kappa),2,1)
#
# all_K <- seq(0,.05,by=.01)/freq
#
# ell <- 0
#
# H <- 8
#
# Pi_t_minus_ell <- exp(ell*model$mu_pi0)
#
# res <- price_Inflation_caps_floors(W = sim.w, # Values of state vector (T x n)
#                                    H = H, # maximum maturity, in model periods
#                                    ell = ell, # Indexation lag
#                                    all_K = all_K, # vector of strikes, expressed at model frequency
#                                    # e.g.: quarter. data and K=.01 means annualized strike of 4%
#                                    Pi_t_minus_ell = Pi_t_minus_ell, # Indexation lag multiplicative factor
#                                    psi = psi.stocks2fact, # Laplace transform of W
#                                    parameterization, # see details below
#                                    max_x = 10000, # settings for Riemann sum comput.
#                                    dx_statio = 10,
#                                    min_dx = 1e-05,
#                                    nb_x1 = 1000)
#
# par(mfrow=c(3,1))
#
# par(plt=c(.1,.95,.1,.8))
#
# plot(w2,type="l",las=1,xlab="",ylab="",col="white",
#      main=expression(paste("(a) Crisis factor (",w[2*','*t],")",sep=""))
# )
# grid()
# lines(w2,lwd=2)
#
# ylim_pi <- seq(-.1,.15,by=.05)
#
# plot(w1,type="l",col="white",xaxt="n",yaxt="n",
#      ylim=c(-.07,.12),
#      las=1,xlab="",ylab="",
#      main=expression(paste("(b) Annualized inflation",sep=""))
# )
# grid()
# axis(side=2,at=ylim_pi,labels = sprintf("%.2f", ylim_pi),las=1,
#      col.axis = "darkgrey",col = "darkgrey")
# lines(annualized_pi_t,col="black",lwd=2)
#
# h <- 2
# plot(res$Caps[,1,h],type="l",las=1,
#      ylim = range(res$Caps[,,h],na.rm = TRUE),
#      col="white",
#      xlab="", ylab="",
#      main=expression(paste("(c) Calls and put prices",sep="")))
# for(i in 1:length(all_K)){
#   lines(res$Caps[,i,h],col=i)
# }
# plot(res$Floors[,1,h],type="l",las=1,
#      ylim = range(res$Floors[,,h],na.rm = TRUE),
#      col="white",
#      xlab="", ylab="",
#      main=expression(paste("(c) Calls and put prices",sep="")))
# for(i in 1:length(all_K)){
#   lines(res$Floors[,i,h],col=i)
# }
#
#
