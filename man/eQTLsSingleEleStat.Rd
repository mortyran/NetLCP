\name{eQTLsSingleEleStat}
\alias{eQTLsSingleEleStat}
\title{single element eQTLs statistics}
\usage{eQTLsSingleEleStat(regData, eQTLsData, filterDegree, selectNode)}
\arguments{
  \item{regData}{the standard output of the function binaryRegulationExtract or multieleRegulation.}
  \item{eQTLsData}{the standard out of the function eQTLsDetection.}
  \item{filterDegree}{an integer for filtering the nodes.}
  \item{selectNode}{select a or a set of nodes to perform statistics.}
}
\value{
 single element eQTLs statistics with bar plot.
}
\description{
single element eQTLs statistics
}
\examples{
 eQTLsSingleEleStat(regData = NULL, eQTLsData = NULL, filterDegree = 100,selectNode = "MIMAT0000423",netLayout = "layout_in_circle")
}
