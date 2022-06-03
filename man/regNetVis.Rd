\name{regNetVis}
\alias{regNetVis}
\title{regulatory Network Visual Function}
\usage{regNetVis(regData, filterDegree, selectCREs, netLayout)}
\arguments{
  \item{regData}{regData standard output of the function binaryRegulation or multieleRegulation.}
  \item{filterDegree}{an integer for filtering the nodes.}
  \item{selectCREs}{a vector of elements or CREs.}
  \item{netLayout}{the layout of network.}
}
\value{
 regulatory network visualization
}
\description{
Regulatory network visualization based on the standard output of the function binaryRegulation or multieleRegulation.
}
\examples{
 regNetVis(regData = NULL, filterDegree = 20, selectCREs = "9978", netLayout = "layout_in_circle")
}
