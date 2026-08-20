\name{DTAM_datasets}
\alias{DTAM_datasets}
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
The datasets are data frames or tibbles. Main contents are:
\describe{
  \item{\code{ACMTermPremium}}{Monthly Adrian-Crump-Moench U.S. Treasury yields,
  term premia, and risk-neutral yields.}
  \item{\code{DAT_GSW}}{Monthly Gurkaynak-Sack-Wright nominal and real
  zero-coupon Treasury yields.}
  \item{\code{DAT_LW}}{Loads the object \code{YC_LW}, selected monthly Liu-Wu
  U.S. zero-coupon Treasury yields.}
  \item{\code{DAT_LW_FULL}}{Loads the object \code{YC_LW_FULL}, the dense monthly
  Liu-Wu U.S. zero-coupon Treasury yield panel from one month through ten years.}
  \item{\code{Data_Macro_EA_quarterly}}{Quarterly euro-area macroeconomic data
  extracted from the Area-Wide Model database.}
  \item{\code{Data_Macro_US_monthly}}{Monthly U.S. macroeconomic variables and
  yield data.}
  \item{\code{Data_Macro_US_quarterly}}{Quarterly U.S. macroeconomic variables
  and yield data.}
  \item{\code{JSTdataset}}{Six-variable, 16-country extract from the
  Jorda-Schularick-Taylor Macrohistory Database used in the book.}
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
source summary with source URLs, bibliographic references, and upstream
licence/terms notes. These notes document the upstream terms known to DTAM; users
remain responsible for checking the current terms of the original data providers.

Main sources and terms:
\describe{
  \item{\code{ACMTermPremium}}{Federal Reserve Bank of New York term-premium
  data. Reference: Adrian, Crump, and Moench (2013), Journal of Financial
  Economics. URL: \url{https://www.newyorkfed.org/research/data_indicators/term-premia-tabs}.
  Reuse is subject to the New York Fed website terms and attribution
  requirements.}
  \item{\code{DAT_GSW}}{Federal Reserve Board nominal and TIPS yield-curve
  estimates. References: Gurkaynak, Sack, and Wright (2007), Journal of Monetary
  Economics; Gurkaynak, Sack, and Wright (2010), American Economic Journal:
  Macroeconomics. URL:
  \url{https://www.federalreserve.gov/data/nominal-yield-curve.htm}. Subject to
  Federal Reserve Board website terms and model disclaimer.}
  \item{\code{DAT_LW}, \code{DAT_LW_FULL}}{Liu-Wu reconstructed Treasury yield
  curves. Reference: Liu and Wu (2021), Journal of Financial Economics. URL:
  \url{https://sites.google.com/view/jingcynthiawu/yield-data}. Cynthia Wu gave
  written permission on 2026-08-19 to include the requested book extract in
  correspondence copied to Yan Liu. Cite the paper. The data are outside
  DTAM's MIT software licence.}
  \item{\code{Data_Macro_EA_quarterly}}{Updated Area-Wide Model Database.
  References: İpek and Kısacıkoğlu (2026), \emph{European
  Economic Review} 181, 105179; Fagan, Henry, and Mestre (2001), ECB Working
  Paper No. 42. URL: \url{https://eabcn.org/data/area-wide-model}. Written
  redistribution permission was received from the database maintainers on
  2026-08-19; they reported EABCN approval. The data are not covered by DTAM's
  MIT software licence.}
  \item{\code{Data_Macro_US_monthly}, \code{Data_Macro_US_quarterly},
  \code{YC_US}, \code{YC_US_weekly}}{Federal Reserve Economic Data (FRED), with
  original source agencies as listed in each FRED series. URL:
  \url{https://fred.stlouisfed.org/}. FRED content is free to access subject to
  FRED terms; third-party series may carry additional restrictions.}
  \item{\code{JSTdataset}}{Six-variable, 16-country extract from the
  Jorda-Schularick-Taylor Macrohistory Database used in the book.
  Reference: Jorda, Schularick, and Taylor (2017), NBER Macroeconomics Annual.
  Users of return data should also cite Jorda et al. (2019), Quarterly Journal of
  Economics; users of bank-balance-sheet ratios should also cite Jorda et al.
  (2021), Review of Economic Studies. URL:
  \url{https://www.macrohistory.net/database/}. The upstream database states a
  Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
  licence; written clarification about bundling this extract is pending.}
  \item{\code{Shiller}}{Robert Shiller online stock-market data, augmented in
  DTAM with FRED series. Reference: Shiller (2015), Irrational Exuberance, 3rd
  edition. URL: \url{https://www.econ.yale.edu/~shiller/data.htm}. Cite the
  source page/book and check the Yale/Shiller page terms before redistribution.}
  \item{\code{SPF}}{Federal Reserve Bank of Philadelphia Survey of Professional
  Forecasters. URL:
  \url{https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/survey-of-professional-forecasters}.
  Reuse is subject to Philadelphia Fed website terms and attribution.}
  \item{\code{YC_Euro}}{ECB euro-area zero-coupon yields from the ECB Data
  Portal. URL: \url{https://data.ecb.europa.eu/}. Publicly released ESCB/ECB
  statistics may be reused free of charge if the source is quoted and the
  statistics/metadata are not modified; third-party data are excluded.}
}
}
\examples{
dtam_datasets()
dtam_datasets(details = TRUE)
YC_US <- dtam_dataset("YC_US")
YC_LW <- dtam_dataset("YC_LW")
}
\keyword{datasets}
