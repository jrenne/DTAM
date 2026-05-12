\name{dtam_dataset}
\alias{dtam_dataset}
\title{Return a DTAM dataset}
\usage{
dtam_dataset(name)
}
\arguments{
\item{name}{Character scalar. Name of the dataset to return. See
\code{\link[=dtam_datasets]{dtam_datasets()}}.}
}
\value{
The requested dataset.
}
\description{
Loads one of the datasets shipped with DTAM and returns it as an object.
}
\examples{
dtam_datasets()
YC_US <- dtam_dataset("YC_US")
}
