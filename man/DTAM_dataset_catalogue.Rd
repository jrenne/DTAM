\name{DTAM_datasets}
\alias{DTAM_datasets}
\alias{ACMTermPremium}
\alias{Credit_spds}
\alias{DAT_GSW}
\alias{DAT_LW}
\alias{DAT_LW_FULL}
\alias{YC_LW}
\alias{YC_LW_FULL}
\alias{Data_Macro_US_monthly}
\alias{Data_Macro_US_quarterly}
\alias{Futures}
\alias{JSTdataset}
\alias{Shiller}
\alias{SPF}
\alias{YC_Euro}
\alias{YC_US}
\alias{YC_US_weekly}
\docType{data}
\title{Datasets shipped with DTAM}
\description{
DTAM includes empirical datasets used in the examples and bookdown code for
discrete-time affine asset-pricing models. The helper function
\code{\link[=dtam_datasets]{dtam_datasets()}} returns a compact catalogue of
these datasets, and \code{\link[=dtam_dataset]{dtam_dataset()}} returns a
dataset by object name. Use
\code{\link[=dtam_dataset_info]{dtam_dataset_info()}} for dimensions, sample
coverage, and variable-level information.
}
\usage{
data(ACMTermPremium)
data(Credit_spds)
data(DAT_GSW)
data(DAT_LW)
data(DAT_LW_FULL)
data(Data_Macro_US_monthly)
data(Data_Macro_US_quarterly)
data(Futures)
data(JSTdataset)
data(Shiller)
data(SPF)
data(YC_Euro)
data(YC_US)
data(YC_US_weekly)
}
\format{
The datasets are data frames or tibbles. Main contents are:
\describe{
  \item{\code{ACMTermPremium}}{Monthly Adrian-Crump-Moench U.S. Treasury yields,
  term premia, and risk-neutral yields.}
  \item{\code{Credit_spds}}{Annual U.S. corporate credit spreads and default
  rates.}
  \item{\code{DAT_GSW}}{Monthly Gurkaynak-Sack-Wright nominal and real
  zero-coupon Treasury yields.}
  \item{\code{DAT_LW}}{Loads the object \code{YC_LW}, selected monthly Liu-Wu
  U.S. zero-coupon Treasury yields.}
  \item{\code{DAT_LW_FULL}}{Loads the object \code{YC_LW_FULL}, the full monthly
  Liu-Wu U.S. zero-coupon Treasury yield panel.}
  \item{\code{Data_Macro_US_monthly}}{Monthly U.S. macroeconomic variables and
  yield data.}
  \item{\code{Data_Macro_US_quarterly}}{Quarterly U.S. macroeconomic variables
  and yield data.}
  \item{\code{Futures}}{Natural gas futures settlement prices.}
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
See \code{\link[=dtam_datasets]{dtam_datasets()}} for a dataset-by-dataset
source summary.
}
\examples{
dtam_datasets()
dtam_datasets(details = TRUE)
YC_US <- dtam_dataset("YC_US")
YC_LW <- dtam_dataset("YC_LW")
}
\keyword{datasets}
