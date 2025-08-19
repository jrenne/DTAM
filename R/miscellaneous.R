
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
