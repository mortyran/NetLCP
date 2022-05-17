\name{eQTLsDetection}
\alias{eQTLsDetection}
\title{Single Biological Elements eQTLs Detection}
\usage{eQTLsDetection(regData)}
\arguments{
  \item{regData}{the standard output from binaryRegulation and multieleRegulation functions,}
}
\value{
The dataframe of eQTLs data with the colname of Elements, SNPID, Location, Source, RegType, DataSource
}
\description{
Detect the eQTLs of single biological elements, which are integrated from ncRNA-eQTLs and eQTLGen.
}
\examples{
eQTLsDetection(regData = NULL)
}
