\name{BioRegElePrioritization}
\alias{BioRegElePrioritization}
\title{Regulatory Elements Prioritization}
\usage{BioRegElePrioritization(transcriptomeList, prioType, empiricalPvalue)}
\arguments{
  \item{transcriptomeList}{a vector consist transcriptome including circRNA, lncRNA, miRNA or mRNA.}
  \item{prioType}{character string naming the prioritizing element types including circRNA, lncRNA, KEGG, Reactome, Wikipathway.}
  \item{empiricalPvalue}{logical value enabling or not enabling empirical p value computing.}
}
\value{
a dataframe of prioritizations.
}
\description{
Prioritize the regulatory elements through decoding the heterogeneous network according to input transcriptome.
}
\examples{
BioRegElePrioritization(transcriptomeList = c("2309", "3838", "MIMAT0000255"),
                                  prioType = "KEGG",
                                  empiricalPvalue = FALSE)
}
