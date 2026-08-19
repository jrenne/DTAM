\name{Data_Macro_EA_quarterly}
\alias{Data_Macro_EA_quarterly}
\docType{data}
\title{Euro-area quarterly macroeconomic data}
\description{
Quarterly euro-area macroeconomic data extracted from the Area-Wide Model
database and formatted for empirical examples in the book.
}
\usage{data(Data_Macro_EA_quarterly)}
\format{
A data frame with the following variables:
\describe{
  \item{date}{Date corresponding to the first month of the quarter.}
  \item{q}{Quarter identifier.}
  \item{YER}{Real GDP.}
  \item{HICP}{Harmonised index of consumer prices.}
  \item{URX}{Unemployment rate.}
  \item{pi}{Quarterly log HICP inflation.}
  \item{dy}{Quarterly log real GDP growth.}
}
}
\source{
Updated Area-Wide Model Database, March 2026 vintage (EA21),
\url{https://eabcn.org/data/area-wide-model}. Cite İpek and Kısacıkoğlu
(2026), \emph{European Economic Review} 181, 105179,
\url{https://doi.org/10.1016/j.euroecorev.2025.105179}.
}
\keyword{datasets}
