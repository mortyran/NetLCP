\name{regVarStat}
\alias{regVarStat}
\title{Regulatory variants Statistics}
\usage{regVarStat(regVar, regulationType, selectNode)}
\arguments{
  \item{regVar}{the standard out of the function regVarDetection.}
  \item{regulationType}{a character representing the regulation data type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".}
  \item{selectNode}{select a or a set of nodes to visulize in the network.}
}
\value{
Regulatory variants Statistics
}
\description{
Regulatory variants Statistics
}
\examples{
regVarStat(regVar = NULL, regulationType = "lncRNA-miRNA-mRNA", selectNode = sampleData)
}
