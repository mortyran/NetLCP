\name{regStat}
\alias{regStat}
\title{biological elements of regulation statistics}
\usage{regStat(regData, filterDegree)}
\arguments{
  \item{regData}{the standard output of the function binaryRegulation or multieleRegulation.}
  \item{filterDegree}{an integer for filtering the nodes.}
  \item{selectNode}{select a or a set of nodes to perform statistics.}
}
\value{
 Elements of regulations statistics with barplot.
}
\description{
biological elements of regulation statistics
}
\examples{
regStat(regData = NULL, filterDegree = 50)
}
