\name{regVarNetVis}
\alias{regVarNetVis}
\title{Variant 'switches' Network Visualization}
\usage{regVarNetVis(regVar, regulationType, selectNode)}
\arguments{
  \item{regVar}{the standard out of the function regVarDetection.}
  \item{regulationType}{a character representing the CREs type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".}
  \item{selectNode}{select a or a set of nodes to visulize in the network.}
}
\value{
Variant 'switches' network visualization
}
\description{
Variant 'switches' network visualization
}
\examples{
regVarNetVis(regVar = NULL, regulationType = "miRNA-mRNA-pathway", selectNode = c("MIMAT0000717","6774", "hsa04550"))
}
