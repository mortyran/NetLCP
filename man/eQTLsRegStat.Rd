\name{eQTLsRegStat}
\alias{eQTLsRegStat}
\title{eQTLs of CREs statistics}
\usage{eQTLsRegStat(regData, eQTLsData, regulationType, topCREs, filterDegree, selectCREs)}
\arguments{
  \item{regData}{the standard output of the function binaryRegulationExtract or multieleRegulation.}
  \item{eQTLsData}{the standard out of the function eQTLsDetection.}
  \item{regulationType}{a character representing the CREs type, it should be one of "circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway".}
  \item{topCREs}{an integer representing the number of the top prioritized CREs to exhibit.}
  \item{filterDegree}{an integer for filtering the nodes.}
  \item{selectCREs}{a vector of elements or CREs.}
}
\value{
 eQTLs on CREs statistics with barplot.
}
\description{
eQTLs of CREs statistics
}
\examples{
 eQTLsRegStat(regData = NULL, eQTLsData = NULL, regulationType = "miRNA-mRNA", topCREs = 10, filterDegree = 40, selectCREs = NULL)
}
