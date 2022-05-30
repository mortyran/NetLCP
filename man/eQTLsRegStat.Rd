\name{eQTLsRegStat}
\alias{eQTLsRegStat}
\title{eQTLs on regulations statistics}
\usage{eQTLsRegStat(regData, eQTLsData, regulationType, filterDegree, selectNode)}
\arguments{
  \item{regData}{the standard output of the function binaryRegulationExtract or multieleRegulation.}
  \item{eQTLsData}{the standard out of the function eQTLsDetection.}
  \item{regulationType}{a character representing the CREs type, it should be one of "circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway".}
  \item{filterDegree}{an integer for filtering the nodes.}
  \item{selectNode}{select a or a set of nodes to perform statistics.}
}
\value{
 eQTLs on regulations statistics with barplot.
}
\description{
eQTLs on regulations statistics
}
\examples{
 eQTLsRegStat(regData = NULL, eQTLsData = NULL, regulationType = "miRNA-mRNA", filterDegree = 40, selectNode = NULL)
}
