
# load libraries:
library(ecb)
library(tidyverse)
library(zoo)
library(fredr)
library(Hmisc)
library(dplyr)


#===============================================================================
# European zero-coupon yields from ECB data
#===============================================================================

maturities <- 1:10
list_of_tickers <- paste("YC.B.U2.EUR.4F.G_N_A.SV_C_YM.SR_",maturities,"Y",sep="")
all_yields <- data.frame()
count <- 0
for(ticker in list_of_tickers){
  count <- count + 1
  x <- get_data(ticker)
  x_reduced <- data.frame(date=x$obstime,
                          yield=x$obsvalue)
  names(x_reduced)[2] <- paste("Y",maturities[count],sep="")
  if(count==1){
    YC_Euro <- x_reduced
  }else{
    YC_Euro <- merge(YC_Euro,x_reduced,by="date")
  }
}
save(YC_Euro,file="data/YC_Euro.rda")



#===============================================================================
# European zero-coupon yields from FRED
#===============================================================================

fredr_set_key("df65e14c054697a52b4511e77fcfa1f3")

f <- function(ticker,freq){
  fredr(series_id = ticker,
        observation_start = start_date,observation_end = end_date,
        frequency = freq,aggregation_method = "avg")
}

list.variables <- c("DTB4WK","DTB3","DTB6",
                    "THREEFY1","THREEFY2","THREEFY3","THREEFY4","THREEFY5",
                    "THREEFY6","THREEFY7","THREEFY8","THREEFY9","THREEFY10")

start_date <- as.Date("1970-01-01")
end_date   <- as.Date("2024-11-01")

freq <- "d"

for(i in 1:length(list.variables)){
  data.var <- f(list.variables[i],freq)
  eval(parse(text = gsub(" ","",paste("data.var.frame = data.frame(date=data.var$date,",
                                      list.variables[i],"=data.var$value)",
                                      sep=""))))
  if(i==1){
    YC_US = data.var.frame
  }else{
    YC_US = merge(YC_US,data.var.frame,by="date",all=TRUE)
  }
}
save(YC_US,file="data/YC_US.rda")



#===============================================================================
# Nominal yield data from GSW 2006
#===============================================================================

download.file("https://www.federalreserve.gov/data/yield-curve-tables/feds200628.csv",
              "data-raw/feds200628.csv")
DAT <- csv.get(file = "data-raw/feds200628.csv", skip = 9)
DAT$Date <- as.Date(DAT$Date)
DAT <- arrange(DAT,Date)
DAT$mdate <- paste(format(DAT$Date,"%Y"),"-",
                   format(DAT$Date,"%m"),sep="")
DAT_mthly <- DAT %>%
  group_by(mdate) %>% #filter(Date == max(Date,na.rm = TRUE))
  summarise_all(function(x){mean(x,na.rm=TRUE)})
DAT_mthly$Date <- as.Date(
  paste(format(DAT_mthly$Date,"%Y"),"-",
        format(DAT_mthly$Date,"%m"),"-01",sep=""))
DAT_GSW_nom <- data.frame(date=DAT_mthly$Date,
                          SVENY01=DAT_mthly$SVENY01,
                          SVENY02=DAT_mthly$SVENY02,
                          SVENY03=DAT_mthly$SVENY03,
                          SVENY04=DAT_mthly$SVENY04,
                          SVENY05=DAT_mthly$SVENY05,
                          SVENY06=DAT_mthly$SVENY06,
                          SVENY07=DAT_mthly$SVENY07,
                          SVENY08=DAT_mthly$SVENY08,
                          SVENY09=DAT_mthly$SVENY09,
                          SVENY10=DAT_mthly$SVENY10,
                          SVENY11=DAT_mthly$SVENY11,
                          SVENY12=DAT_mthly$SVENY12,
                          SVENY13=DAT_mthly$SVENY13,
                          SVENY14=DAT_mthly$SVENY14,
                          SVENY15=DAT_mthly$SVENY15,
                          SVENY16=DAT_mthly$SVENY16,
                          SVENY17=DAT_mthly$SVENY17,
                          SVENY18=DAT_mthly$SVENY18,
                          SVENY19=DAT_mthly$SVENY19,
                          SVENY20=DAT_mthly$SVENY20,
                          SVENY21=DAT_mthly$SVENY21,
                          SVENY22=DAT_mthly$SVENY22,
                          SVENY23=DAT_mthly$SVENY23,
                          SVENY24=DAT_mthly$SVENY24,
                          SVENY25=DAT_mthly$SVENY25,
                          SVENY26=DAT_mthly$SVENY26,
                          SVENY27=DAT_mthly$SVENY27,
                          SVENY28=DAT_mthly$SVENY28,
                          SVENY29=DAT_mthly$SVENY29,
                          SVENY30=DAT_mthly$SVENY30)

save(DAT_GSW_nom,file="data/DAT_GSW_nom.rda")


#===============================================================================
# Real yield data from GSW 2008
#===============================================================================

download.file("https://www.federalreserve.gov/data/yield-curve-tables/feds200805.csv",
              "data-raw/feds200805.csv")
DAT <- csv.get(file = "data-raw/feds200805.csv", skip = 18)
DAT$Date <- as.Date(DAT$Date)
DAT <- arrange(DAT,Date)
DAT$mdate <- paste(format(DAT$Date,"%Y"),"-",
                   format(DAT$Date,"%m"),sep="")
DAT_mthly <- DAT %>%
  group_by(mdate) %>% #filter(Date == max(Date,na.rm = TRUE))
  summarise_all(function(x){mean(x,na.rm=TRUE)})
DAT_mthly$Date <- as.Date(
  paste(format(DAT_mthly$Date,"%Y"),"-",
        format(DAT_mthly$Date,"%m"),"-01",sep=""))
DAT_GSW_real <- data.frame(date=DAT_mthly$Date,
                           TIPSY02=DAT_mthly$TIPSY02,
                           TIPSY03=DAT_mthly$TIPSY03,
                           TIPSY04=DAT_mthly$TIPSY04,
                           TIPSY05=DAT_mthly$TIPSY05,
                           TIPSY06=DAT_mthly$TIPSY06,
                           TIPSY07=DAT_mthly$TIPSY07,
                           TIPSY08=DAT_mthly$TIPSY08,
                           TIPSY09=DAT_mthly$TIPSY09,
                           TIPSY10=DAT_mthly$TIPSY10,
                           TIPSY11=DAT_mthly$TIPSY11,
                           TIPSY12=DAT_mthly$TIPSY12,
                           TIPSY13=DAT_mthly$TIPSY13,
                           TIPSY14=DAT_mthly$TIPSY14,
                           TIPSY15=DAT_mthly$TIPSY15,
                           TIPSY16=DAT_mthly$TIPSY16,
                           TIPSY17=DAT_mthly$TIPSY17,
                           TIPSY18=DAT_mthly$TIPSY18,
                           TIPSY19=DAT_mthly$TIPSY19,
                           TIPSY20=DAT_mthly$TIPSY20)

save(DAT_GSW_real,file="data/DAT_GSW_real.rda")



#===============================================================================
# Data from Cynthia Wu's website
# https://sites.google.com/view/jingcynthiawu/yield-data
# https://docs.google.com/spreadsheets/d/1-wmStGZHLx55dSYi3gQK2vb3F8dMw_Nb/edit?gid=378310471#gid=378310471
#===============================================================================

YC_LW_raw <- csv.get(file = "data-raw/LW_monthly.csv", skip = 8)
maturities_in_month <- matrix(c(1,3,6,9,12*(1:30)),ncol=1)
maturities <- apply(maturities_in_month,1,function(x){ifelse(x<12,
                                                             paste(x,"m",sep=""),
                                                             paste(x/12,"y",sep=""))})
YC_LW <- cbind(YC_LW_raw[,1+maturities_in_month])
names(YC_LW) <- paste("yld_",maturities,sep="")
YC_LW$date <- as.character(YC_LW_raw[,1])
YC_LW$date <- as.Date(paste(substr(YC_LW$date,1,4),"-",substr(YC_LW$date,5,6),"-01",sep=""))

save(YC_LW,file="data/DAT_LW.rda")


#===============================================================================
# Surveys of Professional Forecasters
#===============================================================================

# 1 year CPI -------------------------------------------------------------------
download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_cpi_level.xlsx",
              "data-raw/mean_cpi_level.xlsx")
SPF <- readxl::read_xlsx(path="data-raw/mean_cpi_level.xlsx")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
SPF <- data.frame(date=SPF$date,CPI1=as.numeric(SPF$CPI6))
SPF.CPI <- data.frame(date=SPF$date,CPI1=SPF$CPI1)

# 10 year CPI-------------------------------------------------------------------
download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_cpi10_level.xlsx",
              "data-raw/mean_cpi10_level.xlsx")
SPF <- readxl::read_xlsx(path="data-raw/mean_cpi10_level.xlsx")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
SPF <- data.frame(date=SPF$date,CPI10=as.numeric(SPF$CPI10))
SPF.CPI10<- data.frame(date=SPF$date,CPI10=SPF$CPI10)

# 1 year TBILL -----------------------------------------------------------------
download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_tbill_level.xlsx",
              "data-raw/mean_tbill_level.xlsx")
SPF <- readxl::read_xlsx(path="data-raw/mean_tbill_level.xlsx")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
SPF <- data.frame(date=SPF$date,BILL1=as.numeric(SPF$TBILL6))
SPF.BILL <- data.frame(date=SPF$date,BILL1=SPF$BILL1)

# 10 year TBILL-----------------------------------------------------------------
download.file("https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/mean_bill10_level.xlsx",
              "data-raw/mean_bill10_level.xlsx")
SPF <- readxl::read_xlsx(path="data-raw/mean_bill10_level.xlsx")
SPF$date <- as.Date(paste(SPF$YEAR,"-",1+3*(SPF$QUARTER-1),"-01",sep=""))
SPF <- data.frame(date=SPF$date,BILL10=as.numeric(SPF$BILL10))
SPF.BILL10<- data.frame(date=SPF$date,BILL10=SPF$BILL10)

SPF <- merge(SPF.CPI,SPF.CPI10,by="date",all = TRUE)
SPF <- merge(SPF,SPF.BILL,by="date",all = TRUE)
SPF <- merge(SPF,SPF.BILL10,by="date",all = TRUE)

save(SPF,file="data/SPF.rda")


