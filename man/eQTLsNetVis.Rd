\name{eQTLsNetVis}
\alias{eQTLsNetVis}
\title{eQTLs Regulatory network visualization}
\usage{eQTLsNetVis(regData, eQTLsData, filterDegree, selectCREs, netLayout)}
\arguments{
  \item{regData}{the standard output of the function binaryRegulationExtract or multieleRegulation.}
  \item{eQTLsData}{the standard out of the function eQTLsDetection.}
  \item{filterDegree}{an integer for filtering the nodes.}
  \item{selectCREs}{a vector of elements or CREs.}
  \item{netLayout}{layout of network, including "layout_in_circle" or "layout_nicely".}
}
\value{
 eQTLs Regulatory network visualization
}
\description{
eQTLs Regulatory network visualization based on regData and eQTLsData.
}
\examples{
 eQTLsNetVis(regData = NULL, eQTLsData = NULL, filterDegree = 100, selectCREs = "MIMAT0000423",netLayout = "layout_in_circle")
}
