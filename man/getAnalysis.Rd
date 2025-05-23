% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RFLOMICS-Methods_00_Common.R
\name{getAnalysis}
\alias{getAnalysis}
\alias{getAnalysis,RflomicsMAE-method}
\alias{getAnalysis,RflomicsSE-method}
\alias{getAnalyzedDatasetNames}
\alias{getAnalyzedDatasetNames,RflomicsMAE-method}
\title{Get results from RFLOMICS object}
\usage{
\S4method{getAnalysis}{RflomicsMAE}(object, name = NULL, subName = NULL)

\S4method{getAnalysis}{RflomicsSE}(object, name = NULL, subName = NULL)

\S4method{getAnalyzedDatasetNames}{RflomicsMAE}(object, analyses = NULL)
}
\arguments{
\item{object}{The RflomicsMAE or RflomicsSE object from which to extract
the analysis.}

\item{name}{the name of element to add to metadata slot.}

\item{subName}{the name of sub element to add to metadata slot.}

\item{analyses}{vector of list of analysis name}
}
\value{
The analysis metadata slot (a list of results)
\itemize{
   \item getAnalysis: return list of results from a specific analysis.}
}
\description{
Get a specific analysis results from a RflomicsMAE or a SE.

\itemize{
   \item getAnalyzedDatasetNames: return a list of performed analysis names.}
}
\examples{
# See generateReport for an example that includes getAnalyzedDatasetNames
}
