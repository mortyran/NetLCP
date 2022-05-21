\name{regVarNetVis}
\alias{regVarNetVis}
\title{Regulatory variants Network Visual Function}
\usage{regVarNetVis(regVar, regulationType, selectNode)}
\arguments{
  \item{regVar}{the standard out of the function regVarDetection.}
  \item{regulationType}{a character representing the regulation data type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".}
  \item{selectNode}{select a or a set of nodes to visulize in the network.}
}
\value{
Regulatory variants network visualization
}
\description{
Regulatory variants network visualization
}
\examples{
regVarNetVis(regVar = NULL, regulationType = "miRNA-mRNA-pathway", selectNode = c("MIMAT0000717","6774", "hsa04550"))
}
