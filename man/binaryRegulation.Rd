\name{binaryRegulation}
\alias{binaryRegulation}
\title{Binary CREs in Local Heterogeneous Network}
\usage{binaryRegulation(elmentList, regulationType, allRegulation)}
\arguments{
  \item{elmentList}{a vector consist of elements including circRNA, lncRNA, miRNA ,mRNA or KEGG/Reactome/Wikipathway pathway IDs.}
  \item{regulationType}{character string naming the types of CREs which can be located including "circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway".}
  \item{allRegulation}{a logical value indicates exploring all related CREs from the storage (TRUE) or internal CREs between input elements (FALSE).}
}
\value{
The dataframe of regulatory data with the colnames of Node1, Node2, Datasoure and RegulationType.
}
\description{
Locate the binary CREs in regional Heterogeneous network which are experimentally verified.
}
\examples{
binaryRegulation(elmentList = c("2309", "3838", "MIMAT0000255"), regulationType = "miRNA-mRNA",  allRegulation = FALSE)
}
