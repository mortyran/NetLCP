\name{regVarStat}
\alias{regVarStat}
\title{Variant 'switches' Statistics}
\usage{regVarStat(regVar, regulationType, selectCREs)}
\arguments{
  \item{regVar}{the standard out of the function regVarDetection.}
  \item{regulationType}{a character representing the CREs type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".}
  \item{selectCREs}{a vector of elements or CREs.}
}
\value{
Variant 'switches' Statistics
}
\description{
Variant 'switches' Statistics
}
\examples{
regVarStat(regVar = NULL, regulationType = "lncRNA-miRNA-mRNA", selectCREs = sampleData)
}
