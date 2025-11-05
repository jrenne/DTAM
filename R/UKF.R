

make_sigma_points <- function(w,Sigma,alpha,beta,kappa){
  n <- dim(Sigma)[1]

  lambda <- alpha^2*(n+kappa) - n

  omega <- matrix(0,2*n+1,1)
  omega[1] <- lambda/(lambda + kappa)
  omega[2:(n+1)] <- 0

  Sigma12 <- t(chol(Sigma))

  Chi <- matrix(w,n,2*n+1)
  Chi[,2:(n+1)] <- Chi[,2:(n+1)] + sqrt(lambda + n) * Sigma12
  Chi[,(n+2):(2*n+1)] <- Chi[,(n+2):(2*n+1)] - sqrt(lambda + n) * Sigma12

  return(Chi)
}

# Sigma <- matrix(c(1,.5,.5,1),2,2)
# w <- matrix(0,2,1)
# alpha <- 1
# beta  <- 2
# kappa <- 1
# Chi <- make_sigma_points(w,Sigma,alpha,beta,kappa)
#
# plot(0,0,
#      xlim=range(Chi[1,]),
#      ylim=range(Chi[2,]))
# for(i in 1:dim(Chi)[2]){
#   points(Chi[1,i],Chi[2,i])
# }

