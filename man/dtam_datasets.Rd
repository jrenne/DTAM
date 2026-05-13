\name{dtam_datasets}
\alias{dtam_datasets}
\title{List DTAM datasets}
\usage{
dtam_datasets(details = FALSE)
}
\arguments{
\item{details}{Logical. If \code{TRUE}, add dimensions and date/sample
coverage by loading each dataset.}
}
\value{
A data frame with object names, underlying data names, broad frequency, source,
description, source URL, reference, and upstream licence/terms note. With
\code{details = TRUE}, the data frame also includes the number of observations,
number of variables, and start/end sample markers when a date-like column is
available.
}
\description{
Returns a catalogue of the datasets shipped with DTAM.
}
\examples{
dtam_datasets()
dtam_datasets(details = TRUE)
}
