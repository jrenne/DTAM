#
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
#
#   return(cbind(all_w1,all_w2))
# }
#
# alpha <- .05
# rho   <- .9
# mu    <- 2
# beta  <- rho/mu
# sigma <- .01
# chi   <- .0002
# phi   <- .2
# kappa <- .003
# delta <- .004
#
# model <- list(
#   alpha = alpha,
#   mu    = mu,
#   beta  = beta,
#   sigma = sigma,
#   chi   = chi,
#   phi   = phi,
#   kappa = kappa,
#   delta = delta,
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
# dS <- model$delta + w1 - model$kappa * w2
# S  <- exp(cumsum(dS))
#
# K_over_S <- seq(.8,1.2,by=.2)
#
# parameterization = list()
# parameterization$model <- model
# parameterization$xi0 = 0
# parameterization$xi1 = matrix(0,2,1)
#
# res <- price_Stock_calls_puts(W = sim.w, # Values of state vector (T x n)
#                               S = S, # ex-dividend stock price (T x 1)
#                               H = 5, # maximum maturity, in model periods
#                               a = matrix(c(1,-kappa),ncol=1),
#                               b = delta, # specif. of ex-dividend stock returns
#                               K_over_S = K_over_S, # vector of strikes
#                               psi = psi.stocks2fact, # Laplace transform of W
#                               parameterization = parameterization,
#                               max_x = 1000, # settings for Riemann sum comput.
#                               dx_statio = 1,
#                               min_dx = 1e-05,
#                               nb_x1 = 10000)
#
# par(mfrow=c(3,1))
#
# par(plt=c(.1,.95,.1,.8))
#
# plot(w2,type="l",las=1,xlab="",ylab="",col="white",
#      main=expression(paste("(a) Crisis factor (",w[2*','*t],")",sep=""))
#      )
# grid()
# lines(w2,lwd=2)
#
# plot(w1,type="l",col="white",xaxt="n",yaxt="n",
#      ylim=c(-.15,.15),
#      las=1,xlab="",ylab="",
#      main=expression(paste("(b) Stock returns (in grey)",sep=""))
# )
# grid()
# ylim_dS <- seq(-.15,.15,by=.05)
# axis(side=2,at=ylim_dS,labels = sprintf("%.2f", ylim_dS),las=1,
#      col.axis = "darkgrey",col = "darkgrey")
# lines(dS,col="darkgrey",lwd=2)
#
# par(new=TRUE)
# plot(S,type="l",col="white",xaxt="n",yaxt="n",
#      ylim=c(.2,1.8),
#      xlab="",ylab="")
# lines(S,col="black",lwd=2)
# ylim_S <- seq(0,2,by=.25)
# axis(side=4,at=ylim_S,labels = sprintf("%.2f", ylim_S),las=1,
#      col.axis = "black",col = "black")
#
# plot(res$Puts[,1,5]/S,type="l",las=1,xlab="",ylab="",
#      main=expression(paste("(c) Calls and put prices",sep=""))
# )
#
#
