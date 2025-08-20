
make_recessions <- function(col="#AA55AA44",indic_US=TRUE){
  if(indic_US){
    nber_dates <- nberDates()
    start_recession_dates <- as.Date(as.character(nber_dates[,1]),"%Y%m%d")
    end_recession_dates   <- as.Date(as.character(nber_dates[,2]),"%Y%m%d")
  }else{
    # Euro zone
    start_recession_dates <- as.Date(c("2008-02-01","2011-08-01","2019-11-01"))
    end_recession_dates   <- as.Date(c("2009-05-01","2013-02-01","2020-08-01"))
  }

  for(i in 1:length(start_recession_dates)){
    polygon(c(start_recession_dates[i],start_recession_dates[i],
              end_recession_dates[i],end_recession_dates[i]),
            c(-10000,+10000,+10000,-10000),border=NaN,col=col)
  }
  return(1)
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


autocov <- function(X,n){
  X <- as.matrix(X)
  if(dim(X)[1]<dim(X)[2]){
    X <- t(X)
  }
  T <- dim(X)[1]
  X <- X - matrix(1,T,1) %*% matrix(apply(X,2,mean),nrow=1)
  X.1 <- X[1:(T-n),]
  X.2 <- X[(n+1):T,]
  return(1/T * t(X.1) %*% X.2)
}

NW.LongRunVariance <- function(X,q){
  gamma0 <- autocov(X,0)
  LRV <- gamma0
  weights <- 1 - (1:q)/(q+1)
  for(i in 1:q){
    LRV <- LRV + 2 * weights[i] * autocov(X,i)
  }
  return(LRV)
}
