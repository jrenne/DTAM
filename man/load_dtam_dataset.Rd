\name{load_dtam_dataset}
\alias{load_dtam_dataset}
\title{Load a DTAM dataset into an environment}
\usage{
load_dtam_dataset(name, envir = parent.frame())
}
\arguments{
\item{name}{Character scalar. Name of the dataset to load.}

\item{envir}{Environment in which the dataset should be loaded. Defaults to
the caller environment.}
}
\value{
Invisibly returns the loaded dataset.
}
\description{
Convenience wrapper around \code{\link[utils:data]{data()}} that validates the
dataset name against \code{\link[=dtam_datasets]{dtam_datasets()}}.
}
\examples{
load_dtam_dataset("YC_Euro")
}
