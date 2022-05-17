\name{regVarDetection}
\alias{regVarDetection}
\title{Mutations on regulations Detection}
\usage{regVarDetection(regVar, regulationType)}
\arguments{
  \item{regData}{the dataframe of output from binaryRegulation and multieleRegulation functions,
                 accepting miRNA-mRNA, miRNA-mRNA-pathway, lncRNA-miRNA-mRNA, circRNA-miRNA-mRNA, lncRNA-miRNA-mRNA-pathway, circRNA-miRNA-mRNA-pathway.}
  \item{regulationType}{a character representing the extracted regulation data type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".}
}
\value{
Regulatory Variants Detection.
}
\description{
Detect the Mutations on regulations, which are integrated from LnCeVar and miRNASNPv3
}
\examples{
regVarDetection(regVar = NULL, regulationType = "lncRNA-miRNA-mRNA", selectNode = sampleData)
}
