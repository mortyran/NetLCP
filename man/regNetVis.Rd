\name{regNetVis}
\alias{regNetVis}
\title{regulatory Network Visual Function}
\usage{regNetVis(regData, filterDegree, selectNode, netLayout)}
\arguments{
  \item{regData}{regData standard output of the function binaryRegulation or multieleRegulation.}
  \item{filterDegree}{an integer for filtering the nodes.}
  \item{selectNode}{select a or a set of nodes to visulize in the network.}
  \item{netLayout}{the layout of network.}
}
\value{
 regulatory network visualization
}
\description{
Regulatory network visualization based on the standard output of the function binaryRegulation or multieleRegulation.
}
\examples{
 regNetVis(regData = NULL, filterDegree = 20, selectNode = "9978", netLayout = "layout_in_circle")
}
