\name{regVarStat}
\alias{regVarStat}
\title{Variant 'switches' Statistics}
\usage{regVarStat(regVar, regulationType, selectNode)}
\arguments{
  \item{regVar}{the standard out of the function regVarDetection.}
  \item{regulationType}{a character representing the CREs type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".}
  \item{selectNode}{select a or a set of nodes to visulize in the network.}
}
\value{
Variant 'switches' Statistics
}
\description{
Variant 'switches' Statistics
}
\examples{
regVarStat(regVar = NULL, regulationType = "lncRNA-miRNA-mRNA", selectNode = sampleData)
}
