\name{DTAM_data_objects}
\alias{ACMTermPremium}
\alias{DAT_GSW}
\alias{DAT_LW}
\alias{DAT_LW_FULL}
\alias{YC_LW}
\alias{YC_LW_FULL}
\alias{Data_Macro_US_monthly}
\alias{Data_Macro_US_quarterly}
\alias{JSTdataset}
\alias{Shiller}
\alias{SPF}
\alias{YC_Euro}
\alias{YC_US}
\alias{YC_US_weekly}
\docType{data}
\title{DTAM data objects}
\description{
Empirical data objects shipped with DTAM. Use
\code{\link[=dtam_datasets]{dtam_datasets()}} to list available datasets and
\code{\link[=dtam_dataset_info]{dtam_dataset_info()}} to retrieve dimensions,
sample coverage, variables, source URLs, references, and upstream licence/terms
notes.
}
\usage{
data(ACMTermPremium)
data(DAT_GSW)
data(DAT_LW)
data(DAT_LW_FULL)
data(Data_Macro_EA_quarterly)
data(Data_Macro_US_monthly)
data(Data_Macro_US_quarterly)
data(JSTdataset)
data(Shiller)
data(SPF)
data(YC_Euro)
data(YC_US)
data(YC_US_weekly)
}
\format{
The objects are data frames or tibbles:
\describe{
  \item{\code{ACMTermPremium}}{Monthly Adrian-Crump-Moench U.S. Treasury yields,
  term premia, and risk-neutral yields.}
  \item{\code{DAT_GSW}}{Monthly Gurkaynak-Sack-Wright nominal and real
  zero-coupon Treasury yields.}
  \item{\code{DAT_LW}}{Loads the object \code{YC_LW}, selected monthly Liu-Wu
  U.S. zero-coupon Treasury yields.}
  \item{\code{DAT_LW_FULL}}{Loads the object \code{YC_LW_FULL}, the full monthly
  Liu-Wu U.S. zero-coupon Treasury yield panel.}
  \item{\code{Data_Macro_EA_quarterly}}{Quarterly euro-area macroeconomic data
  extracted from the Area-Wide Model database.}
  \item{\code{Data_Macro_US_monthly}}{Monthly U.S. macroeconomic variables and
  yield data.}
  \item{\code{Data_Macro_US_quarterly}}{Quarterly U.S. macroeconomic variables
  and yield data.}
  \item{\code{JSTdataset}}{Jorda-Schularick-Taylor Macrohistory Database.}
  \item{\code{Shiller}}{Monthly Shiller stock-market, valuation, price, and
  consumption data.}
  \item{\code{SPF}}{Philadelphia Fed Survey of Professional Forecasters data
  for CPI, GDP, and Treasury-bill forecasts.}
  \item{\code{YC_Euro}}{Daily euro-area zero-coupon yields from the ECB.}
  \item{\code{YC_US}}{Daily U.S. Treasury yields from FRED.}
  \item{\code{YC_US_weekly}}{Weekly U.S. rates, policy targets, and Treasury
  yields from FRED.}
}
}
\source{
Run \code{dtam_datasets()} for compact source metadata, or
\code{dtam_dataset_info("YC_LW")} for one dataset. A longer source and licence
summary is installed at \code{system.file("DATA_SOURCES.md", package = "DTAM")}.
}
\examples{
YC_LW <- dtam_dataset("YC_LW")
dtam_dataset_info("YC_LW")
}
\keyword{datasets}
