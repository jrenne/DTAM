#' Add recession shading to a time-series plot
#'
#' Draws shaded polygons corresponding to U.S. or euro-area recession dates on
#' the current base-R plot.
#'
#' @param col Fill color used for the shaded areas.
#' @param indic_US Logical. If `TRUE`, use U.S. recession dates from `tseries`;
#'   otherwise use the hard-coded euro-area dates included in the package.
#' @param MIN,MAX Vertical bounds used for the polygon.
#'
#' @return Invisibly returns `NULL`.
#'
#' @examples
#' plot(as.Date(c("2007-01-01", "2012-01-01")),
#'      c(0, 1), type = "n", xlab = "", ylab = "")
#' make_recessions(col = "#CCCCCC66", MIN = 0, MAX = 1)
#'
#' @export
make_recessions <- function(col="#AA55AA44",indic_US=TRUE,
                            MIN=-10000,MAX=+10000){
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
            c(MIN,MAX,MAX,MIN),border=NaN,col=col)
  }
  invisible(NULL)
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
#' Fast Gaussian multi-horizon aggregation
#'
#' Computes the mean, autoregressive loading, and covariance associated with the
#' cumulated Gaussian VAR state `x_{t+1} + ... + x_{t+h}` for a collection of
#' horizons.
#'
#' @param model List containing `mu`, `Phi`, and `Sigma` for the Gaussian
#'   VAR(1).
#' @param Maturities_decompo Integer decomposition matrix describing the target
#'   maturities.
#'
#' @return A list with arrays `MU`, `PHI`, and `SIGMA`, plus a vector `H` of
#'   maturities.
#'
#' @examples
#' model <- list(
#'   mu = matrix(c(0, 0.1), 2, 1),
#'   Phi = matrix(c(0.8, 0.1,
#'                  0.0, 0.7), 2, 2, byrow = TRUE),
#'   Sigma = diag(c(0.2, 0.1))
#' )
#' decomp <- matrix(c(2, 0,
#'                    1, 2), 2, 2, byrow = TRUE)
#' res <- compute_Gaussian_fast(model, decomp)
#' res$H
#'
#' @export
