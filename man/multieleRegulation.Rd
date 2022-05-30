\name{multieleRegulation}
\alias{multieleRegulation}
\title{Multielements CREs in Local Heterogeneous Network}
\usage{multieleRegulation(elementList, regulationType, allRegulation)}
\arguments{
  \item{elementList}{a vector consist of elements including circRNA, lncRNA, miRNA, mRNA or pathway IDs.}
  \item{regulationType}{character string naming the type of CREs including "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA", "lncRNA-miRNA-pathway" and "circRNA-miRNA-pathway".}
  \item{allRegulation}{a logical value indicates locating all related CREs from the storage or internal CREs between input elements.}
}
\value{
The dataframe of regulatory data with the colnames of Node1, Node2, Datasoure and RegulationType.
}
\description{
Locate the multielements CREs in regional heterogeneous network which are experimentally verified.
}
\examples{
multieleRegulation(elementList = c("2309", "3838", "ENSG00000259366", "MIMAT0000255"), regulationType = "lncRNA-miRNA-mRNA",  allRegulation = FALSE)
}
