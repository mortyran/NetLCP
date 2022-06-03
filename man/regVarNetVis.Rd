\name{regVarNetVis}
\alias{regVarNetVis}
\title{Variant 'switches' Network Visualization}
\usage{regVarNetVis(regVar, regulationType, selectCREs)}
\arguments{
  \item{regVar}{the standard out of the function regVarDetection.}
  \item{regulationType}{a character representing the CREs type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".}
  \item{selectCREs}{a vector of elements or CREs.}
}
\value{
Variant 'switches' network visualization
}
\description{
Variant 'switches' network visualization
}
\examples{
regVarNetVis(regVar = NULL, regulationType = "miRNA-mRNA-pathway", selectCREs = c("MIMAT0000717","6774", "hsa04550"))
}
