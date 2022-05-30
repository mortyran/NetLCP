\name{regVarDetection}
\alias{regVarDetection}
\title{Variant 'switches' in CREs Detection}
\usage{regVarDetection(regVar, regulationType)}
\arguments{
  \item{regData}{the dataframe of output from binaryRegulation and multieleRegulation functions,
                 accepting miRNA-mRNA, miRNA-mRNA-pathway, lncRNA-miRNA-mRNA, circRNA-miRNA-mRNA, lncRNA-miRNA-mRNA-pathway, circRNA-miRNA-mRNA-pathway.}
  \item{regulationType}{a character representing the CREs type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".}
}
\value{
Variant 'switches' Detection.
}
\description{
Detect the Variant 'switches' in CREs.
}
\examples{
regVarDetection(regVar = NULL, regulationType = "lncRNA-miRNA-mRNA", selectNode = sampleData)
}
