\name{dtam_dataset_info}
\alias{dtam_dataset_info}
\title{Describe a DTAM dataset}
\usage{
dtam_dataset_info(name)
}
\arguments{
\item{name}{Character scalar. Name of the dataset to describe. See
\code{\link[=dtam_datasets]{dtam_datasets()}}.}
}
\value{
A list with two elements: \code{metadata}, a one-row data frame with
dataset-level information, and \code{variables}, a data frame with variable
names, classes, and missing-value counts.
}
\description{
Returns metadata, dimensions, sample coverage, and variable-level information
for one dataset shipped with DTAM.
}
\examples{
info <- dtam_dataset_info("Data_Macro_EA_quarterly")
info$metadata
info$variables
}
