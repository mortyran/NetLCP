\name{eQTLsNetVis}
\alias{eQTLsNetVis}
\title{eQTLs Regulatory network visualization}
\usage{eQTLsNetVis(regData, eQTLsData, filterDegree, selectNode, netLayout)}
\arguments{
  \item{regData}{the standard output of the function binaryRegulationExtract or multieleRegulation.}
  \item{eQTLsData}{the standard out of the function eQTLsDetection.}
  \item{filterDegree}{an integer for filtering the nodes.}
  \item{selectNode}{select a or a set of nodes to visulize in the network.}
  \item{netLayout}{the layout of network.}
}
\value{
 eQTLs Regulatory network visualization
}
\description{
eQTLs Regulatory network visualization based on regData and eQTLsData.
}
\examples{
 eQTLsNetVis(regData = NULL, eQTLsData = NULL, filterDegree = 100,selectNode = "MIMAT0000423",netLayout = "layout_in_circle")
}
