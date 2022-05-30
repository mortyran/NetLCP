\name{regStat}
\alias{regStat}
\title{biological elements of CREs statistics}
\usage{regStat(regData, filterDegree)}
\arguments{
  \item{regData}{the standard output of the function binaryRegulation or multieleRegulation.}
  \item{filterDegree}{an integer for filtering the nodes.}
  \item{selectNode}{select a or a set of nodes to perform statistics.}
}
\value{
 Elements of CREs statistics with barplot.
}
\description{
biological elements of CREs statistics
}
\examples{
regStat(regData = NULL, filterDegree = 50)
}
